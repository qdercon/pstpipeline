// -----------------------------------------------------------------------------
// Q-learning model for PST training data + affect + overall/block time
//------------------------------------------------------------------------------
// References:
//// https://www.pnas.org/doi/10.1073/pnas.1407535111
//// https://osf.io/9g2zw (replication of the above with Stan code)
//// https://psyarxiv.com/v86bx (Beta distribution response)
//// https://psyarxiv.com/bwv58 (passage-of-time-dysphoria, inspired time pars)
//------------------------------------------------------------------------------

data {
  int<lower=1> N, T;          // # participants, max # of trials
  array[N] int Tsubj;         // # of trials for acquisition phase

  array[N, T] int option1;    // LHS option (1-6)
  array[N, T] int option2;    // RHS option (1-6)
  array[N, T] int choice;     // choice (1 = chose option 1, 3, or 5)
  matrix[N, T] reward;        // coded as 1 (reward) or -1 (no reward)

  matrix[N, T] affect;        // includes 0 and 1, needs to be transformed
  array[N, T] int question;   // from 1 to 3 (happy, confident, engaged)
  matrix[N, T] ovl_time;      // in hours to keep weights relatively small
  matrix[N, T] blk_time;      // in hours to keep weights relatively small
}

transformed data {
  // Default values to initialize the vector of expected values
  vector[6] initial_values;
  initial_values = rep_vector(0, 6);

  int s;
  matrix[N, T] affect_tr;
  s = N*T;
  for (t in 1:T) {
    affect_tr[:, t] = ((affect[:, t]*(s-1))+0.5)/s; // Smith & Verkuilen, 2006
  }
}

parameters {
  // group-level RL parameters
  vector[2] mu_ql;
  vector<lower=0>[2] sigma_ql;

  // group-level weights
  matrix[5, 3] mu_wt;
  matrix<lower=0>[5, 3] sigma_wt;

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
  matrix[N, 3] w1_b_pr;
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
  matrix[N, 3] w1_b;
  matrix[N, 3] w2;
  matrix[N, 3] w3;
  matrix<lower=0, upper=1>[N, 3] gamma;
  matrix[N, 3] phi;

  for (p in 1:3) {
    w0[:, p]    = mu_wt[1, p] + sigma_wt[1, p] * w0_pr[:, p];
    w1_o[:, p]  = mu_wt[2, p] + sigma_wt[2, p] * w1_o_pr[:, p];
    w1_b[:, p]  = mu_wt[3, p] + sigma_wt[3, p] * w1_b_pr[:, p];
    w2[:, p]    = mu_wt[4, p] + sigma_wt[4, p] * w2_pr[:, p];
    w3[:, p]    = mu_wt[5, p] + sigma_wt[5, p] * w3_pr[:, p];
    gamma[:, p] = Phi_approx(mu_gm[p] + sigma_gm[p] * gm_pr[:, p]);
    phi[:, p]   = exp(aff_mu_phi[p] + aff_sigma_phi[p] * phi_pr[:, p]);
  }
}

model {
  // hyperpriors on QL parameters
  mu_ql    ~ normal(0, 1);
  sigma_ql ~ normal(0, 0.2);

  // hyperpriors on the weights
  for (p in 1:3) {
    mu_wt[:, p]    ~ normal(0, 1);
    sigma_wt[:, p] ~ exponential(0.1);
  }

  // hyperpriors on gamma
  mu_gm    ~ normal(0, 1);
  sigma_gm ~ exponential(0.1);

  // hyperpriors on the beta distribution precision
  aff_mu_phi    ~ normal(0, 1);
  aff_sigma_phi ~ exponential(0.1);

  // priors on QL parameters
  alpha_pr ~ normal(0, 1);
  beta_pr  ~ normal(0, 1);

  // priors on the weights + gamma + beta distribution precision
  for (p in 1:3) {
   w0_pr[:, p]   ~ normal(0, 1);
   w1_o_pr[:, p] ~ normal(0, 1);
   w1_b_pr[:, p] ~ normal(0, 1);
   w2_pr[:, p]   ~ normal(0, 1);
   w3_pr[:, p]   ~ normal(0, 1);
   gm_pr[:, p]   ~ normal(0, 1);
   phi_pr[:, p]  ~ normal(0, 1);
  }

  for (i in 1:N) {
    int co;                       // Chosen option
    real pe;                      // Prediction error
    real delta;                   // Difference between two options
    vector[6] ev;                 // Expected values per symbol
    vector[Tsubj[i]] decay_vec;   // Weighting of previous trials

    real aff_mu_cond;             // Conditional mean of the beta distribution
    vector[Tsubj[i]] shape_a;     // Beta distribution shape parameter alpha
    vector[Tsubj[i]] shape_b;     // Beta distribution shape parameter beta

    row_vector[Tsubj[i]] ev_vec;
    row_vector[Tsubj[i]] pe_vec;

    ev     = initial_values;
    ev_vec = rep_row_vector(0, Tsubj[i]);
    pe_vec = rep_row_vector(0, Tsubj[i]);

    // Acquisition Phase
    for (t in 1:Tsubj[i]) {

      co = (choice[i, t] > 0) ? option1[i, t] : option2[i, t];

      // Luce choice rule
      delta = ev[option1[i, t]] - ev[option2[i, t]];
      choice[i, t] ~ bernoulli_logit(beta[i] * delta);

      pe = reward[i, t] - ev[co];
      pe_vec[t] = pe;

      ev[co] += alpha[i] * pe;
      ev_vec[t] = ev[co];

      decay_vec[t] = pow(gamma[i, question[i, t]], t - 1);

      aff_mu_cond = inv_logit(
        w0[i, question[i, t]] +
        w1_o[i, question[i, t]] * ovl_time[i, t] +
        w1_b[i, question[i, t]] * blk_time[i, t] +
        w2[i, question[i, t]] * (reverse(ev_vec[:t]) * decay_vec[:t]) +
        w3[i, question[i, t]] * (reverse(pe_vec[:t]) * decay_vec[:t])
      );

      // add machine precision to ensure shape parameters > 0
      shape_a[t] =
        aff_mu_cond * phi[i, question[i, t]] + machine_precision();
      shape_b[t] =
        phi[i, question[i, t]] * (1-aff_mu_cond) + machine_precision();
    }

    // increment log density (vectorised)
    affect_tr[i, :Tsubj[i]] ~ beta(shape_a, shape_b);
  }
}

generated quantities {
  // For group-level parameters
  real<lower=0,upper=1>  mu_alpha;
  real<lower=0,upper=10> mu_beta;

  vector[3] mu_w0;
  vector[3] mu_w1_o;
  vector[3] mu_w1_b;
  vector[3] mu_w2;
  vector[3] mu_w3;
  vector[3] mu_gamma;

  // For log-likelihood calculation
  vector[N] log_lik;
  matrix[N, T] y_pred; // will get NaNs for those not predicted

  for (i in 1:N) {
    for (t in 1:Tsubj[i]) {
      y_pred[i, t] = -1;
    }
  };

  mu_alpha = Phi_approx(mu_ql[1]);
  mu_beta  = Phi_approx(mu_ql[2]) * 10;

  for (p in 1:3) {
    mu_w0[p]    = mu_wt[1, p];
    mu_w1_o[p]  = mu_wt[2, p];
    mu_w1_b[p]  = mu_wt[3, p];
    mu_w2[p]    = mu_wt[4, p];
    mu_w3[p]    = mu_wt[5, p];
    mu_gamma[p] = Phi_approx(mu_gm[p]);
  }

  for (i in 1:N) {
    int co;                       // Chosen option
    real pe;                      // Prediction error
    real delta;                   // Difference between two options
    vector[6] ev;                 // Expected values per symbol
    vector[Tsubj[i]] decay_vec;   // Weighting of previous trials

    real aff_mu_cond;             // Conditional mean of the beta distribution
    vector[Tsubj[i]] shape_a;     // Beta distribution shape parameter alpha
    vector[Tsubj[i]] shape_b;     // Beta distribution shape parameter beta

    row_vector[Tsubj[i]] ev_vec;
    row_vector[Tsubj[i]] pe_vec;

    ev     = initial_values;
    ev_vec = rep_row_vector(0, Tsubj[i]);
    pe_vec = rep_row_vector(0, Tsubj[i]);

    log_lik[i] = 0;

    // Acquisition Phase
    for (t in 1:Tsubj[i]) {
      co = (choice[i, t] > 0) ? option1[i, t] : option2[i, t];

      // Luce choice rule
      delta = ev[option1[i, t]] - ev[option2[i, t]];
      log_lik[i] += bernoulli_logit_lpmf(choice[i, t] | beta[i] * delta);

      pe = reward[i, t] - ev[co];
      pe_vec[t] = pe;

      ev[co] += alpha[i] * pe;
      ev_vec[t] = ev[co];

      decay_vec[t] = pow(gamma[i, question[i, t]], t - 1);

      aff_mu_cond = inv_logit(
        w0[i, question[i, t]] +
        w1_o[i, question[i, t]] * ovl_time[i, t] +
        w1_b[i, question[i, t]] * blk_time[i, t] +
        w2[i, question[i, t]] * (reverse(ev_vec[:t]) * decay_vec[:t]) +
        w3[i, question[i, t]] * (reverse(pe_vec[:t]) * decay_vec[:t])
      );

      // add machine precision to ensure shape parameters > 0
      shape_a[t] =
        aff_mu_cond * phi[i, question[i, t]] + machine_precision();
      shape_b[t] =
        phi[i, question[i, t]] * (1-aff_mu_cond) + machine_precision();

      // generate posterior predictions
      y_pred[i, t] = beta_rng(shape_a[t], shape_b[t]);
    }
     // increment log likelihood with affect model
    log_lik[i] += beta_lpdf(affect_tr[i, :Tsubj[i]] | shape_a, shape_b);
  }
}
