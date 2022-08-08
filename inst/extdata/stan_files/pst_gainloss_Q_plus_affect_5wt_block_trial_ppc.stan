// -----------------------------------------------------------------------------
// Gain-loss Q-learning model for PST training data with affect data extension
//------------------------------------------------------------------------------
// References:
//// https://www.pnas.org/doi/10.1073/pnas.1407535111
//// https://osf.io/9g2zw (replication of the above with Stan code)
//// https://psyarxiv.com/bwv58 (passage-of-time-dysphoria, inspired time pars)
//------------------------------------------------------------------------------

data {
  int<lower=1> N, T;          // # participants, max # of trials
  int<lower=1> Tsubj[N];      // # of trials for acquisition phase

  int<lower=-1,upper=6> option1[N, T];
  int<lower=-1,upper=6> option2[N, T];
  int<lower=-1,upper=1> choice[N, T];
  real reward[N, T];

  real<lower=0, upper=1> affect[N, T];
  int<lower=0, upper=3> question[N, T];
  real<lower=0> ovl_trial[N, T];
  real<lower=0> blk_trial[N, T];
}

transformed data {
  // Default values to initialize the vector of expected values
  vector[6] initial_values;
  initial_values = rep_vector(0, 6);
}

parameters {
  // Group-level parameters
  vector[3] mu_ql_pr;
  vector<lower=0>[3] sigma_ql;

  real<lower=2> nu; // minimum 2 df

  // weights
  real mu_wt[5, 3];
  real<lower=0> sigma_wt[5, 3];

  /* decay factor (gamma)
  real<lower=0> alpha_g[3];
  real<lower=0> beta_g[3];
  */

  // scale parameter (sigma_t)
  real<lower=0> alpha_s;
  real<lower=0> beta_s;

  // individual-level parameters for Matt trick
  vector[N] alpha_pos_pr;
  vector[N] alpha_neg_pr;
  vector[N] beta_pr;

  // individual-level weights
  real w0[N, 3];
  real w1_o[N, 3];
  real w1_b[N, 3];
  real w2[N, 3];
  real w3[N, 3];

  real<lower=0, upper=1> gamma[N, 3]; // forgetting factor
  vector<lower=0>[N] sigma_t; // scale of t distribution
}

transformed parameters {
  vector<lower=0,upper=1>[N] alpha_pos;
  vector<lower=0,upper=1>[N] alpha_neg;
  vector<lower=0,upper=10>[N] beta;

  alpha_pos = Phi_approx(mu_ql_pr[1] + sigma_ql[1] * alpha_pos_pr);
  alpha_neg = Phi_approx(mu_ql_pr[2] + sigma_ql[2] * alpha_neg_pr);
  beta      = Phi_approx(mu_ql_pr[3] + sigma_ql[3] * beta_pr) * 10;
}

model {
  // hyperpriors on QL parameters
  mu_ql_pr ~ normal(0, 1);
  sigma_ql ~ normal(0, 0.2);

  // priors on QL parameters
  alpha_pos_pr ~ normal(0, 1);
  alpha_neg_pr ~ normal(0, 1);
  beta_pr      ~ normal(0, 1);

  // priors on the weights
  for (p in 1:3) {
    w0[:, p]   ~ normal(mu_wt[1, p], sigma_wt[1, p]);
    w1_o[:, p] ~ normal(mu_wt[2, p], sigma_wt[2, p]);
    w1_b[:, p] ~ normal(mu_wt[3, p], sigma_wt[3, p]);
    w2[:, p]   ~ normal(mu_wt[4, p], sigma_wt[4, p]);
    w3[:, p]   ~ normal(mu_wt[5, p], sigma_wt[5, p]);
  }

  // priors on the t distribution
  sigma_t  ~ gamma(alpha_s, beta_s);
  nu       ~ exponential(0.1);

  /*
  for (p in 1:3) {
    gamma[:, p] ~ beta(alpha_g[p], beta_g[p]);
  }*/

  // priors on gamma - not including hyperpriors improves fit.
  for (p in 1:3) {
    gamma[:, p] ~ beta(1, 1);
  }

  // hyperpriors on the weights
  for (p in 1:3) {
    mu_wt[:, p]    ~ normal(0, 1);
    sigma_wt[:, p] ~ uniform(0.001, 1000);
  }

  // hyperpriors on the t distribution
  alpha_s  ~ uniform(0.001, 1000);
  beta_s   ~ uniform(0.001, 1000);

  /*// hyperpriors on gamma
  alpha_g  ~ gamma(1, 0.1);
  beta_g   ~ gamma(1, 0.1);*/

  for (i in 1:N) {
    int co;                       // Chosen option
    real delta;                   // Difference between two options
    real pe;                      // Prediction error
    real alpha;
    vector[6] ev;                 // Expected values per symbol
    vector[Tsubj[i]] decayvec;
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

      decayvec[t] = pow(gamma[i, question[i, t]], t - 1);

      affect[i, t] ~ student_t(
        nu,
        w0[i, question[i, t]] +
        w1_o[i, question[i, t]] * ovl_trial[i, t] +
        w1_b[i, question[i, t]] * blk_trial[i, t] +
        w2[i, question[i, t]] * (reverse(ev_vec[:t]) * decayvec[:t]) +
        w3[i, question[i, t]] * (reverse(pe_vec[:t]) * decayvec[:t]),
        sigma_t[i]
      );
    }
  }
}

generated quantities {
  // For group-level parameters
  real<lower=0,upper=1>  mu_alpha_pos;
  real<lower=0,upper=1>  mu_alpha_neg;
  real<lower=0,upper=10> mu_beta;

  real mu_w0[3];
  real mu_w1_o[3];
  real mu_w1_b[3];
  real mu_w2[3];
  real mu_w3[3];
  real mu_q_gamma[3];

  // For log-likelihood calculation
  real log_lik[N];
  real y_pred[N, T]; // will get NaNs for those not predicted

  for (i in 1:N) {
    for (t in 1:Tsubj[i]) {
      y_pred[i, t] = -1;
    }
  };

  mu_alpha_pos = Phi_approx(mu_ql_pr[1]);
  mu_alpha_neg = Phi_approx(mu_ql_pr[2]);
  mu_beta      = Phi_approx(mu_ql_pr[3]) * 10;

  for (p in 1:3) {
    mu_w0[p]    = Phi_approx(mu_wt[1, p]);
    mu_w1_o[p]  = Phi_approx(mu_wt[2, p]);
    mu_w1_b[p]  = Phi_approx(mu_wt[3, p]);
    mu_w2[p]    = Phi_approx(mu_wt[4, p]);
    mu_w3[p]    = Phi_approx(mu_wt[5, p]);
    mu_q_gamma[p] = quantile(gamma[p], 0.5);
  }

  {
    for (i in 1:N) {
      int co;                       // Chosen option
      int aff_num;
      real delta;                   // Difference between two options
      real pe;                      // Prediction error
      real alpha;
      vector[6] ev;                 // Expected values per symbol
      vector[Tsubj[i]] decayvec;
      row_vector[Tsubj[i]] ev_vec;
      row_vector[Tsubj[i]] pe_vec;

      ev         = initial_values;
      ev_vec     = rep_row_vector(0, Tsubj[i]);
      pe_vec     = rep_row_vector(0, Tsubj[i]);
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

        decayvec[t] = pow(gamma[i, question[i, t]], t - 1);

        log_lik[i] += student_t_lpdf(affect[i, t] |
          nu,
          w0[i, question[i, t]] +
          w1_o[i, question[i, t]] * ovl_trial[i, t] +
          w1_b[i, question[i, t]] * blk_trial[i, t] +
          w2[i, question[i, t]] * (reverse(ev_vec[:t]) * decayvec[:t]) +
          w3[i, question[i, t]] * (reverse(pe_vec[:t]) * decayvec[:t]),
          sigma_t[i]
        );

        y_pred[i, t] = student_t_rng(
          nu,
          w0[i, question[i, t]] +
          w1_o[i, question[i, t]] * ovl_trial[i, t] +
          w1_b[i, question[i, t]] * blk_trial[i, t] +
          w2[i, question[i, t]] * (reverse(ev_vec[:t]) * decayvec[:t]) +
          w3[i, question[i, t]] * (reverse(pe_vec[:t]) * decayvec[:t]),
          sigma_t[i]
        );
      }
    }
  }
}
