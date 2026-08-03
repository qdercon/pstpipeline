// -----------------------------------------------------------------------------
// Gain-loss Q-learning model for PST training data + change in affect rating
//------------------------------------------------------------------------------
// References:
//// https://www.pnas.org/doi/10.1073/pnas.1407535111
//// https://osf.io/9g2zw (replication of the above with Stan code)
//// https://www.nature.com/articles/s41562-023-01519-7 (mood drift over time)
//------------------------------------------------------------------------------

functions {
  /* Choice regressor plus per-trial (pre-update) EV and PE for one
     participant. Storing ev BEFORE this trial's alpha * pe update avoids the
     partial collinearity with pe that the lag terms in affect_eta rely on. */
  matrix pst_regressors(int ti,
                        array[] int option1_i,
                        array[] int option2_i,
                        array[] int choice_i,
                        row_vector reward_i,
                        real alpha_pos_i,
                        real alpha_neg_i) {
    matrix[ti, 3] out; // [ , 1] = delta, [ , 2] = ev_vec, [ , 3] = pe_vec
    vector[6] ev = rep_vector(0, 6);

    for (t in 1:ti) {
      int co = (choice_i[t] > 0) ? option1_i[t] : option2_i[t];
      real pe;
      real a;

      // Luce choice rule (i.e., EVs of non-seen options assumed not to matter)
      out[t, 1] = ev[option1_i[t]] - ev[option2_i[t]];

      pe = reward_i[t] - ev[co];

      out[t, 2] = ev[co]; // pre-update EV
      out[t, 3] = pe;

      a = (pe >= 0) ? alpha_pos_i : alpha_neg_i;
      ev[co] += a * pe;
    }
    return out;
  }

  /* Linear predictor for the affect model: baseline + overall-time slope,
     plus a lag-weighted sum of EV/PE over the trials intervening since the
     previous rating of the same question (including the current trial). */
  vector affect_eta(int ti,
                    array[] int question_i,
                    array[] int int_trials_i,
                    row_vector ovl_time_i,
                    vector ev_vec,
                    vector pe_vec,
                    row_vector w0_i,
                    row_vector w1_o_i,
                    matrix w2_i,
                    matrix w3_i) {
    vector[ti] eta;
    for (t in 1:ti) {
      real e = w0_i[question_i[t]] + w1_o_i[question_i[t]] * ovl_time_i[t];
      for (j in 1:int_trials_i[t]) {
        int t1 = t + 1 - j;
        e += w2_i[question_i[t1], j] * ev_vec[t1]
           + w3_i[question_i[t1], j] * pe_vec[t1];
      }
      eta[t] = e;
    }
    return eta;
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
                        array[,] int int_trials,
                        array[] row_vector affect_tr,
                        array[] row_vector ovl_time,
                        vector alpha_pos,
                        vector alpha_neg,
                        vector beta,
                        matrix w0,
                        matrix w1_o,
                        array[] matrix w2,
                        array[] matrix w3,
                        matrix phi) {
    real lp = 0;
    for (n in 1:size(seq_i)) {
      int i  = seq_i[n];
      int ti = Tsubj[i];
      array[ti] int qn = question[i, 1:ti];

      matrix[ti, 3] reg = pst_regressors(ti, option1[i], option2[i], choice[i],
                                         reward[i], alpha_pos[i], alpha_neg[i]);

      vector[ti] eta = affect_eta(ti, qn, int_trials[i, 1:ti],
                                  ovl_time[i][1:ti], reg[ : , 2], reg[ : , 3],
                                  w0[i], w1_o[i], w2[i], w3[i]);

      vector[ti] mu = inv_logit(eta);
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
  int<lower=1> N, T, I;              // # participants, max # of trials, max (theoretical) # of intervening trials
  array[N] int<lower=1,upper=T> Tsubj;  // # of trials for acquisition phase

  array[N, T] int<lower=1,upper=6> option1;  // LHS option (1-6)
  array[N, T] int<lower=1,upper=6> option2;  // RHS option (1-6)
  array[N, T] int<lower=0,upper=1> choice;   // choice (1 = chose option 1)
  matrix[N, T] reward;                       // coded as 1 (reward) or -1 (no)

  array[N] row_vector<lower=0,upper=1>[T] affect;  // includes 0 and 1
  array[N, T] int<lower=1,upper=3> question;       // 1 happy, 2 conf, 3 engag
  array[N] row_vector<lower=0>[T] ovl_time;        // in hours
  array[N, T] int<lower=0,upper=I> int_trials;     // intervening trials between ratings
}

transformed data {
  array[N] int seq_i;
  for (i in 1:N) seq_i[i] = i;

  // transform affect to be strictly between 0 and 1 (Smithson & Verkuilen, 2006)
  array[N] row_vector[T] affect_tr;
  for (i in 1:N) {
    affect_tr[i] = ((affect[i] * (N - 1)) + 0.5) / N;
  }
}

parameters {
  // group-level RL parameters
  vector[3] mu_ql;
  vector<lower=0>[3] sigma_ql;

  // group-level weights
  matrix[3, 4] mu_wt; // 3 questions x 4 weights
  matrix<lower=0>[3, 4] sigma_wt;

  // group-level beta distribution precision (phi)
  vector[3] aff_mu_phi;
  vector<lower=0>[3] aff_sigma_phi;

  // individual-level RL parameters
  vector[N] alpha_pos_pr;
  vector[N] alpha_neg_pr;
  vector[N] beta_pr;

  // individual-level weights + forgetting factor
  matrix[N, 3] w0_pr;
  matrix[N, 3] w1_o_pr;

  // time-dependent weights (2nd level, by individual)
  matrix[N, 3] w2_i_pr;
  matrix[N, 3] w3_i_pr;
  matrix<lower=0>[N, 3] sigma_w2_pr;
  matrix<lower=0>[N, 3] sigma_w3_pr;

  // individual-level time-point weight parameters (3rd level)
  array[N] matrix[3, I] w2_pr;
  array[N] matrix[3, I] w3_pr;

  // individual-level beta distribution precision parameter
  matrix[N, 3] phi_pr;
}

transformed parameters {
  vector<lower=0, upper=1>[N] alpha_pos;
  vector<lower=0, upper=1>[N] alpha_neg;
  vector<lower=0, upper=10>[N] beta;

  alpha_pos = Phi_approx(mu_ql[1] + sigma_ql[1] * alpha_pos_pr);
  alpha_neg = Phi_approx(mu_ql[2] + sigma_ql[2] * alpha_neg_pr);
  beta      = Phi_approx(mu_ql[3] + sigma_ql[3] * beta_pr) * 10;

  matrix[N, 3] w0;
  matrix[N, 3] w1_o;
  matrix[N, 3] w2_i;
  matrix[N, 3] w3_i;
  matrix[N, 3] phi;

  for (q in 1:3) {
    w0[:, q]   = mu_wt[q, 1] + sigma_wt[q, 1] * w0_pr[:, q];
    w1_o[:, q] = mu_wt[q, 2] + sigma_wt[q, 2] * w1_o_pr[:, q];
    w2_i[:, q] = mu_wt[q, 3] + sigma_wt[q, 3] * w2_i_pr[:, q];
    w3_i[:, q] = mu_wt[q, 4] + sigma_wt[q, 4] * w3_i_pr[:, q];
    phi[:, q]  = exp(aff_mu_phi[q] + aff_sigma_phi[q] * phi_pr[:, q]);
  }

  // time-dependent weights
  array[N] matrix[3, I] w2;
  array[N] matrix[3, I] w3;

  // get individuals' weights for each trial back
  for (i in 1:N) {
    for (q in 1:3) { // over all lags I
      w2[i, q, :] = w2_i[i, q] + sigma_w2_pr[i, q] * w2_pr[i, q, :];
      w3[i, q, :] = w3_i[i, q] + sigma_w3_pr[i, q] * w3_pr[i, q, :];
    }
  }
}

model {
  // hyperpriors on QL parameters
  mu_ql    ~ normal(0, 1);
  sigma_ql ~ exponential(1);

  // hyperpriors on the weights
  for (q in 1:3) {
    mu_wt[q]     ~ normal(0, 1);
    sigma_wt[q]  ~ exponential(1);
  }

  // hyperpriors on the beta distribution precision
  aff_mu_phi    ~ normal(0, 1);
  aff_sigma_phi ~ exponential(1);

  // priors on QL parameters
  alpha_pos_pr ~ normal(0, 2);
  alpha_neg_pr ~ normal(0, 2);
  beta_pr      ~ normal(0, 1);

  // priors on the mean-level weights and beta distribution precision
  for (q in 1:3) {
    w0_pr[:, q]       ~ normal(0, 1);
    w1_o_pr[:, q]     ~ normal(0, 1);
    w2_i_pr[:, q]     ~ normal(0, 1);
    w3_i_pr[:, q]     ~ normal(0, 1);
    sigma_w2_pr[:, q] ~ exponential(1);
    sigma_w3_pr[:, q] ~ exponential(1);
    phi_pr[:, q]      ~ normal(0, 1);
    for (j in 1:I) {
      w2_pr[:, q, j]  ~ normal(0, 1);
      w3_pr[:, q, j]  ~ normal(0, 1);
    }
  }

  if (!prior_only) {
    target += reduce_sum(partial_sum_lupmf, seq_i, grainsize,
                         Tsubj, option1, option2, choice, reward, question,
                         int_trials, affect_tr, ovl_time,
                         alpha_pos, alpha_neg, beta, w0, w1_o, w2, w3, phi);
  }
}

generated quantities {
  // group-level parameter means
  real<lower=0,upper=1>  mu_alpha_pos = Phi_approx(mu_ql[1]);
  real<lower=0,upper=1>  mu_alpha_neg = Phi_approx(mu_ql[2]);
  real<lower=0,upper=10> mu_beta      = Phi_approx(mu_ql[3]) * 10;

  vector[3] mu_w0  = mu_wt[ : , 1];
  vector[3] mu_w1_o = mu_wt[ : , 2];
  vector[3] mu_phi  = exp(aff_mu_phi);

  array[3] int bsl  = { 1, 1, 2 };
  array[3] int comp = { 2, 3, 3 };

  vector[3] phi_diff = mu_phi[bsl] - mu_phi[comp];

  // posterior predictions and log-likelihoods
  vector[run_gq ? N : 0] log_lik;
  array[run_gq ? N : 0] row_vector[T] y_pred;

  if (run_gq) {
    for (i in 1:N) {
      real ll_choice_i = 0;
      real ll_affect_i = 0;
      int ti = Tsubj[i];
      array[ti] int qn = question[i, 1:ti];

      matrix[ti, 3] reg = pst_regressors(ti, option1[i], option2[i], choice[i],
                                         reward[i], alpha_pos[i], alpha_neg[i]);

      vector[ti] eta = affect_eta(ti, qn, int_trials[i, 1:ti],
                                  ovl_time[i][1:ti], reg[ : , 2], reg[ : , 3],
                                  w0[i], w1_o[i], w2[i], w3[i]);

      vector[ti] mu = inv_logit(eta);
      vector[ti] ph = to_vector(phi[i, qn]);
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
