// -----------------------------------------------------------------------------
// Gain-loss Q-learning model for PST training data with affect data extension
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
  vector[N] grp1;              // whether or not participant was non-distanced
  vector[N] grp2;              // whether or not participant was distanced

  array[N,T] int option1;
  array[N,T] int option2;
  array[N,T] int choice;
  matrix[N,T] reward;

  matrix[N, T] affect;        // includes 0 and 1, needs to be transformed
  array[N, T] int question;   // from 1 to 3 (happy, confident, engaged)
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
  // group-level RL parameters (for each group)
  vector[3] mu_ql_g1;
  vector<lower=0>[3] sigma_ql_g1;
  vector[3] mu_ql_g2;
  vector<lower=0>[3] sigma_ql_g2;

  // group-level weights
  matrix[3, 3] mu_wt_g1;
  matrix<lower=0>[3, 3] sigma_wt_g1;
  matrix[3, 3] mu_wt_g2;
  matrix<lower=0>[3, 3] sigma_wt_g2;

  // group-level parameters for decay factor (gamma)
  vector[3] mu_gm_g1;
  vector<lower=0>[3] sigma_gm_g1;
  vector[3] mu_gm_g2;
  vector<lower=0>[3] sigma_gm_g2;

  // group-level beta distribution precision (phi)
  vector[3] mu_phi_g1;
  vector<lower=0>[3] sigma_phi_g1;
  vector[3] mu_phi_g2;
  vector<lower=0>[3] sigma_phi_g2;

  // individual-level RL parameters
  vector[N] alpha_pos_pr;
  vector[N] alpha_neg_pr;
  vector[N] beta_pr;

  // individual-level weights + forgetting factor
  matrix[N, 3] w0_pr;
  matrix[N, 3] w2_pr;
  matrix[N, 3] w3_pr;
  matrix[N, 3] gm_pr;

  // individual-level affect precision (phi)
  matrix[N, 3] phi_pr;

}

transformed parameters {
  vector<lower=0, upper=1>[N] alpha_pos;
  vector<lower=0, upper=1>[N] alpha_neg;
  vector<lower=0, upper=10>[N] beta;

  alpha_pos = Phi_approx(
    (mu_ql_g1[1] * grp1 + sigma_ql_g1[1] * (alpha_pos_pr .* grp1)) +
    (mu_ql_g2[1] * grp2 + sigma_ql_g2[1] * (alpha_pos_pr .* grp2))
  );
  alpha_neg = Phi_approx(
    (mu_ql_g1[2] * grp1 + sigma_ql_g1[2] * (alpha_neg_pr .* grp1)) +
    (mu_ql_g2[2] * grp2 + sigma_ql_g2[2] * (alpha_neg_pr .* grp2))
  );
  beta = Phi_approx(
    (mu_ql_g1[3] * grp1 + sigma_ql_g1[3] * (beta_pr .* grp1)) +
    (mu_ql_g2[3] * grp2 + sigma_ql_g2[3] * (beta_pr .* grp2))
  ) * 10;

  matrix[N, 3] w0;
  matrix[N, 3] w2;
  matrix[N, 3] w3;
  matrix<lower=0, upper=1>[N, 3] gamma;
  matrix[N, 3] phi;

  for (p in 1:3) {
    w0[:, p] =
      (mu_wt_g1[1, p] * grp1 + sigma_wt_g1[1, p] * (w0_pr[:, p] .* grp1)) +
      (mu_wt_g2[1, p] * grp2 + sigma_wt_g2[1, p] * (w0_pr[:, p] .* grp2));
    w2[:, p] =
      (mu_wt_g1[2, p] * grp1 + sigma_wt_g1[2, p] * (w2_pr[:, p] .* grp1)) +
      (mu_wt_g2[2, p] * grp2 + sigma_wt_g2[2, p] * (w2_pr[:, p] .* grp2));
    w3[:, p] =
      (mu_wt_g1[3, p] * grp1 + sigma_wt_g1[3, p] * (w3_pr[:, p] .* grp1)) +
      (mu_wt_g2[3, p] * grp2 + sigma_wt_g2[3, p] * (w3_pr[:, p] .* grp2));

    gamma[:, p] = Phi_approx(
      (mu_gm_g1[p] * grp1 + sigma_gm_g1[p] * (gm_pr[:, p] .* grp1)) +
      (mu_gm_g2[p] * grp2 + sigma_gm_g2[p] * (gm_pr[:, p] .* grp2))
    );
    phi[:, p] = exp(
      (mu_phi_g1[p] * grp1 + sigma_phi_g1[p] * (phi_pr[:, p] .* grp1)) +
      (mu_phi_g2[p] * grp2 + sigma_phi_g2[p] * (phi_pr[:, p] .* grp2))
    );
  }
}

model {
  // hyperpriors on QL parameters
  mu_ql_g1    ~ normal(0, 1);
  mu_ql_g2    ~ normal(0, 1);
  sigma_ql_g1 ~ normal(0, 0.2);
  sigma_ql_g2 ~ normal(0, 0.2);

  // hyperpriors on the weights
  for (p in 1:3) {
    mu_wt_g1[:, p]    ~ normal(0, 1);
    mu_wt_g2[:, p]    ~ normal(0, 1);
    sigma_wt_g1[:, p] ~ exponential(0.1);
    sigma_wt_g2[:, p] ~ exponential(0.1);
  }

  // hyperpriors on gamma
  mu_gm_g1     ~ normal(0, 1);
  mu_gm_g2     ~ normal(0, 1);
  sigma_gm_g1  ~ exponential(0.1);
  sigma_gm_g2  ~ exponential(0.1);

  // hyperpriors on the beta distribution precision
  mu_phi_g1    ~ normal(0, 1);
  mu_phi_g2    ~ normal(0, 1);
  sigma_phi_g1 ~ exponential(0.1);
  sigma_phi_g2 ~ exponential(0.1);

  // priors on QL parameters
  alpha_pos_pr ~ normal(0, 1);
  alpha_neg_pr ~ normal(0, 1);
  beta_pr      ~ normal(0, 1);

  // priors on the weights + gamma + beta distribution precision
  for (p in 1:3) {
   w0_pr[:, p]  ~ normal(0, 1);
   w2_pr[:, p]  ~ normal(0, 1);
   w3_pr[:, p]  ~ normal(0, 1);
   gm_pr[:, p]  ~ normal(0, 1);
   phi_pr[:, p] ~ normal(0, 1);
  }

  for (i in 1:N) {
    int co;                       // Chosen option
    real pe;                      // Prediction error
    real alpha;                   // Learning rate (positive or negative)
    real delta;                   // Difference between two options
    vector[6] ev;                 // Expected values per symbol
    vector[Tsubj[i]] decay_vec;   // Weighting of previous trials

    real mu_cond;             // Conditional mean of the beta distribution
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

      alpha = (pe >= 0) ? alpha_pos[i] : alpha_neg[i];
      ev[co] += alpha * pe;
      ev_vec[t] = ev[co];

      decay_vec[t] = pow(gamma[i, question[i, t]], t - 1);

      mu_cond = inv_logit(
        w0[i, question[i, t]] +
        w2[i, question[i, t]] * (reverse(ev_vec[:t]) * decay_vec[:t]) +
        w3[i, question[i, t]] * (reverse(pe_vec[:t]) * decay_vec[:t])
      );

      // add machine precision to ensure shape parameters > 0
      shape_a[t] =
        mu_cond * phi[i, question[i, t]] + machine_precision();
      shape_b[t] =
        phi[i, question[i, t]] * (1-mu_cond) + machine_precision();
    }

    // increment log density (vectorised)
    affect_tr[i, :Tsubj[i]] ~ beta(shape_a, shape_b);
  }
}

generated quantities {
  // For group-level parameters
  vector[3] mu_alpha_pos;
  vector[3] mu_alpha_neg;
  vector[3] mu_beta;

  matrix[3, 3] mu_w0;
  matrix[3, 3] mu_w2;
  matrix[3, 3] mu_w3;
  matrix[3, 3] mu_gamma;

  // For log-likelihood calculation
  vector[N] log_lik;
  matrix[N, T] y_pred; // will get NaNs for those not predicted

  for (i in 1:N) {
    for (t in 1:Tsubj[i]) {
      y_pred[i, t] = -1;
    }
  };

  mu_alpha_pos[1] = Phi_approx(mu_ql_g1[1]);
  mu_alpha_pos[2] = Phi_approx(mu_ql_g2[1]);
  mu_alpha_pos[3] = Phi_approx(mu_ql_g2[1] - mu_ql_g1[1]);

  mu_alpha_neg[1] = Phi_approx(mu_ql_g1[2]);
  mu_alpha_neg[2] = Phi_approx(mu_ql_g2[2]);
  mu_alpha_neg[3] = Phi_approx(mu_ql_g2[2] - mu_ql_g1[2]);

  mu_beta[1]      = Phi_approx(mu_ql_g1[3]) * 10;
  mu_beta[2]      = Phi_approx(mu_ql_g2[3]) * 10;
  mu_beta[3]      = Phi_approx(mu_ql_g2[3] - mu_ql_g1[3]) * 10;

  for (p in 1:3) {
    mu_w0[1, p] = mu_wt_g1[1, p];
    mu_w0[2, p] = mu_wt_g2[1, p];
    mu_w0[3, p] = mu_wt_g2[1, p] - mu_wt_g1[1, p];

    mu_w2[1, p] = mu_wt_g1[2, p];
    mu_w2[2, p] = mu_wt_g2[2, p];
    mu_w2[3, p] = mu_wt_g2[2, p] - mu_wt_g1[2, p];

    mu_w3[1, p] = mu_wt_g1[3, p];
    mu_w3[2, p] = mu_wt_g2[3, p];
    mu_w3[3, p] = mu_wt_g2[3, p] - mu_wt_g1[3, p];

    mu_gamma[1, p] = Phi_approx(mu_gm_g1[p]);
    mu_gamma[2, p] = Phi_approx(mu_gm_g2[p]);
    mu_gamma[3, p] = Phi_approx(mu_gm_g2[p] - mu_gm_g1[p]);
  }

  for (i in 1:N) {
    int co;                       // Chosen option
    real pe;                      // Prediction error
    real alpha;                   // Learning rate (positive or negative)
    real delta;                   // Difference between two options
    vector[6] ev;                 // Expected values per symbol
    vector[Tsubj[i]] decay_vec;   // Weighting of previous trials

    real mu_cond;             // Conditional mean of the beta distribution
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
      alpha = (pe >= 0) ? alpha_pos[i] : alpha_neg[i];

      ev[co] += alpha * pe;
      ev_vec[t] = ev[co];

      decay_vec[t] = pow(gamma[i, question[i, t]], t - 1);

      mu_cond = inv_logit(
        w0[i, question[i, t]] +
        w2[i, question[i, t]] * (reverse(ev_vec[:t]) * decay_vec[:t]) +
        w3[i, question[i, t]] * (reverse(pe_vec[:t]) * decay_vec[:t])
      );

      // add machine precision to ensure shape parameters > 0
      shape_a[t] =
        mu_cond * phi[i, question[i, t]] + machine_precision();
      shape_b[t] =
        phi[i, question[i, t]] * (1-mu_cond) + machine_precision();

      // generate posterior predictions
      y_pred[i, t] = beta_rng(shape_a[t], shape_b[t]);
    }
    // increment log likelihood with affect model
    log_lik[i] += beta_lpdf(affect_tr[i, :Tsubj[i]] | shape_a, shape_b);
  }
}
