// -----------------------------------------------------------------------------
// Q-learning model for PST training data + affect + block number
//------------------------------------------------------------------------------
// References:
//// https://www.pnas.org/doi/10.1073/pnas.1407535111
//// https://osf.io/9g2zw (replication of the above with Stan code)
//// https://psyarxiv.com/v86bx (Beta distribution response)
//// https://psyarxiv.com/bwv58 (passage-of-time-dysphoria, inspired time pars)
//------------------------------------------------------------------------------

functions {
  matrix pst_regressors(int ti,
                        array[] int option1_i,
                        array[] int option2_i,
                        array[] int choice_i,
                        row_vector reward_i,
                        array[] int question_i,
                        real alpha_i,
                        row_vector gamma_i) {
    matrix[ti, 3] out;
    vector[6] ev        = rep_vector(0, 6);      // expected values per symbol
    row_vector[3] ev_sum = rep_row_vector(0, 3); // summed EVs by question
    row_vector[3] pe_sum = rep_row_vector(0, 3); // summed PEs by question

    for (t in 1:ti) {
      int co = (choice_i[t] > 0) ? option1_i[t] : option2_i[t];
      int qn = question_i[t];
      real pe;

      // Luce choice rule (i.e., EVs of non-seen options assumed not to matter)
      out[t, 1] = ev[option1_i[t]] - ev[option2_i[t]];

      pe = reward_i[t] - ev[co];

      ev_sum += ev[co];
      pe_sum += pe;

      ev[co] += alpha_i * pe;

      // store summed EVs and PEs for this trial
      out[t, 2] = ev_sum[qn];
      out[t, 3] = pe_sum[qn];

      // decay EVs and PEs (i.e., gamma weighted sum over prev. trials)
      ev_sum = ev_sum .* gamma_i;
      pe_sum = pe_sum .* gamma_i;
    }
    return out;
  }

  /* Conditional mean of the Beta distribution for one participant. */
  vector affect_mu(matrix reg,
                   array[] int question_i,
                   row_vector w0_i,
                   row_vector w1_o_i,
                   row_vector w2_i,
                   row_vector w3_i,
                   vector block_no_i) {
    return inv_logit(to_vector(w0_i[question_i])
                     + to_vector(w1_o_i[question_i]) .* block_no_i
                     + to_vector(w2_i[question_i])   .* reg[ : , 2]
                     + to_vector(w3_i[question_i])   .* reg[ : , 3]);
  }

  /* Beta shape parameters from the conditional mean mu and precision phi
     (Smithson & Verkuilen, 2006). Offset by machine precision so both shapes
     stay strictly positive as mu -> 0/1 or phi -> 0; beta_proportion_lpdf
     cannot guarantee this internally, which is what causes the frequent
     rejected proposals under ADVI/HMC. */
  matrix affect_shapes(vector mu_i, vector phi_i) {
    int ti = num_elements(mu_i);
    matrix[ti, 2] out;
    out[ : , 1] = mu_i .* phi_i + machine_precision();
    out[ : , 2] = (1 - mu_i) .* phi_i + machine_precision();
    return out;
  }

  /* Sliced over participants for reduce_sum. */
  real partial_sum_lpmf(array[] int seq_i, int start, int end,
                        array[] int Tsubj,
                        array[,] int option1,
                        array[,] int option2,
                        array[,] int choice,
                        matrix reward,
                        array[,] int question,
                        array[] row_vector affect_tr,
                        array[] vector block_no_tr,
                        vector alpha,
                        vector beta,
                        matrix w0,
                        matrix w1_o,
                        matrix w2,
                        matrix w3,
                        matrix gamma,
                        matrix phi) {
    real lp = 0;
    for (n in 1:size(seq_i)) {
      int i  = seq_i[n];
      int ti = Tsubj[i];
      array[ti] int qn = question[i, 1:ti];

      matrix[ti, 3] reg = pst_regressors(ti, option1[i], option2[i], choice[i],
                                         reward[i], question[i],
                                         alpha[i], gamma[i]);

      vector[ti] mu = affect_mu(reg, qn, w0[i], w1_o[i], w2[i], w3[i],
                                block_no_tr[i][1:ti]);
      matrix[ti, 2] shapes = affect_shapes(mu, to_vector(phi[i, qn]));

      lp += bernoulli_logit_lupmf(choice[i, 1:ti] | beta[i] * reg[ : , 1]);
      lp += beta_lupdf(affect_tr[i][1:ti] | shapes[ : , 1], shapes[ : , 2]);
    }
    return lp;
  }
}

data {
  int<lower=0,upper=1> run_gq;      // 1 to run generated quantities, 0 to skip
  int<lower=0,upper=1> prior_only;  // 1 to sample the prior (for prior PCs)
  int<lower=1> grainsize;           // reduce_sum grainsize; 1 = let Stan choose
  int<lower=1> N, T;                // # participants, max # of trials
  array[N] int<lower=1,upper=T> Tsubj;  // # of trials for acquisition phase

  array[N, T] int<lower=1,upper=6> option1;  // LHS option (1-6)
  array[N, T] int<lower=1,upper=6> option2;  // RHS option (1-6)
  array[N, T] int<lower=0,upper=1> choice;   // choice (1 = chose option 1)
  matrix[N, T] reward;                       // coded as 1 (reward) or -1 (no)

  array[N] row_vector<lower=0,upper=1>[T] affect;  // includes 0 and 1
  array[N, T] int<lower=1,upper=3> question;       // 1 happy, 2 conf, 3 engag
  array[N, T] int<lower=1,upper=6> block_no;       // from 1 to 6
}

transformed data {
  array[N] int seq_i;
  for (i in 1:N) seq_i[i] = i;

  // transform affect to be strictly between 0 and 1 (Smithson & Verkuilen, 2006)
  array[N] row_vector[T] affect_tr;
  for (i in 1:N) {
    affect_tr[i] = ((affect[i] * (N - 1)) + 0.5) / N;
  }

  // transpose block number for efficient indexing
  array[N] vector[T] block_no_tr;
  for (i in 1:N) {
    block_no_tr[i] = to_vector(block_no[i]);
  }
}

parameters {
  // group-level RL parameters
  vector[2] mu_ql;
  vector<lower=0>[2] sigma_ql;

  // group-level weights
  matrix[3, 4] mu_wt; // 3 questions x 4 weights
  matrix<lower=0>[3, 4] sigma_wt;

  // group-level parameters for decay factor (gamma)
  vector[3] mu_gm;
  vector<lower=0>[3] sigma_gm;

  // group-level beta distribution precision (phi)
  vector[3] aff_mu_phi;
  vector<lower=0>[3] aff_sigma_phi;

  // individual-level RL parameters
  vector[N] alpha_pr;
  vector[N] beta_pr;

  // individual-level weights + forgetting factor
  matrix[N, 3] w0_pr;
  matrix[N, 3] w1_o_pr;
  matrix[N, 3] w2_pr;
  matrix[N, 3] w3_pr;
  matrix[N, 3] gm_pr;

  // individual-level affect precision (phi)
  matrix[N, 3] phi_pr;
}

transformed parameters {
  vector<lower=0,upper=1>[N] alpha;
  vector<lower=0,upper=10>[N] beta;

  alpha = Phi_approx(mu_ql[1] + sigma_ql[1] * alpha_pr);
  beta  = Phi_approx(mu_ql[2] + sigma_ql[2] * beta_pr) * 10;

  matrix[N, 3] w0;
  matrix[N, 3] w1_o;
  matrix[N, 3] w2;
  matrix[N, 3] w3;
  matrix<lower=0, upper=1>[N, 3] gamma;
  matrix[N, 3] phi;

  for (q in 1:3) {
    w0[:, q]    = mu_wt[q, 1] + sigma_wt[q, 1] * w0_pr[:, q];
    w1_o[:, q]  = mu_wt[q, 2] + sigma_wt[q, 2] * w1_o_pr[:, q];
    w2[:, q]    = mu_wt[q, 3] + sigma_wt[q, 3] * w2_pr[:, q];
    w3[:, q]    = mu_wt[q, 4] + sigma_wt[q, 4] * w3_pr[:, q];
    gamma[:, q] = Phi_approx(mu_gm[q] + sigma_gm[q] * gm_pr[:, q]);
    phi[:, q]   = exp(aff_mu_phi[q] + aff_sigma_phi[q] * phi_pr[:, q]);
  }
}

model {
  // hyperpriors on QL parameters
  mu_ql    ~ normal(0, 1);
  sigma_ql ~ normal(0, 0.2);

  // hyperpriors on the weights
  for (q in 1:3) {
    mu_wt[q]    ~ normal(0, 1);
    sigma_wt[q] ~ exponential(1);
  }

  // hyperpriors on gamma
  mu_gm    ~ normal(0, 1);
  sigma_gm ~ exponential(1);

  // hyperpriors on the beta distribution precision
  aff_mu_phi    ~ normal(0, 1);
  aff_sigma_phi ~ exponential(1);

  // priors on QL parameters
  alpha_pr ~ normal(0, 1);
  beta_pr  ~ normal(0, 1);

  // priors on the weights + gamma + beta distribution precision
  for (q in 1:3) {
    w0_pr[:, q]   ~ normal(0, 1);
    w1_o_pr[:, q] ~ normal(0, 1);
    w2_pr[:, q]   ~ normal(0, 1);
    w3_pr[:, q]   ~ normal(0, 1);
    gm_pr[:, q]   ~ normal(0, 1);
    phi_pr[:, q]  ~ normal(0, 1);
  }

  if (!prior_only) {
    target += reduce_sum(partial_sum_lupmf, seq_i, grainsize,
                         Tsubj, option1, option2, choice, reward, question,
                         affect_tr, block_no_tr, alpha, beta,
                         w0, w1_o, w2, w3, gamma, phi);
  }
}

generated quantities {
  // group-level parameter means
  real<lower=0,upper=1>  mu_alpha = Phi_approx(mu_ql[1]);
  real<lower=0,upper=10> mu_beta  = Phi_approx(mu_ql[2]) * 10;

  vector[3] mu_w0    = mu_wt[ : , 1];
  vector[3] mu_w1_o  = mu_wt[ : , 2];
  vector[3] mu_w2    = mu_wt[ : , 3];
  vector[3] mu_w3    = mu_wt[ : , 4];
  vector[3] mu_gamma = Phi_approx(mu_gm);
  vector[3] mu_phi   = exp(aff_mu_phi);

  // difference in weights between questions
  array[3] int bsl  = { 1, 1, 2 };
  array[3] int comp = { 2, 3, 3 };

  vector[3] w0_diff    = mu_w0[bsl]    - mu_w0[comp];
  vector[3] w1_o_diff  = mu_w1_o[bsl]  - mu_w1_o[comp];
  vector[3] w2_diff    = mu_w2[bsl]    - mu_w2[comp];
  vector[3] w3_diff    = mu_w3[bsl]    - mu_w3[comp];
  vector[3] gamma_diff = mu_gamma[bsl] - mu_gamma[comp];
  vector[3] phi_diff   = mu_phi[bsl]   - mu_phi[comp];

  // posterior predictions and log-likelihoods
  vector[run_gq ? N : 0] log_lik;
  array[run_gq ? N : 0] row_vector[T] y_pred;

  // calculate log-likelihoods and posterior predictions
  if (run_gq) {
    for (i in 1:N) {
      real ll_choice_i = 0;
      real ll_affect_i = 0;
      int ti = Tsubj[i];
      array[ti] int qn = question[i, 1:ti];

      matrix[ti, 3] reg = pst_regressors(ti, option1[i], option2[i], choice[i],
                                         reward[i], question[i],
                                         alpha[i], gamma[i]);

      vector[ti] mu     = affect_mu(reg, qn, w0[i], w1_o[i], w2[i], w3[i],
                                    block_no_tr[i][1:ti]);
      vector[ti] ph     = to_vector(phi[i, qn]);
      matrix[ti, 2] shapes = affect_shapes(mu, ph);

      ll_choice_i = bernoulli_logit_lpmf(choice[i, 1:ti] | beta[i] * reg[ : , 1]);
      ll_affect_i = beta_lpdf(affect_tr[i][1:ti] | shapes[ : , 1], shapes[ : , 2]);
      log_lik[i] = ll_choice_i + ll_affect_i;

      // posterior predictions, back-transformed to the original scale
      y_pred[i] = rep_row_vector(-1, T);
      y_pred[i][1:ti] =
        ((to_row_vector(beta_rng(shapes[ : , 1], shapes[ : , 2])) * N) - 0.5)
        / (N - 1);
    }
  }
}
