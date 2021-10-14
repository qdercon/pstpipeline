// Gain-loss Q-learning model for PST test data (posterior predictive checks)
data {
  int<lower=1> N;             // Number of subjects
  int<lower=1> T;             // Maximum # of trials
  int<lower=1> T_t;           // Maximum # of test trials
  int<lower=1> Tsubj[N];      // # of trials for acquisition phase
  int<lower=1> Tsubj_t[N];    // # of trials for test phase

  int<lower=-1,upper=6> option1[N, T];
  int<lower=-1,upper=6> option2[N, T];
  int<lower=-1,upper=1> choice[N, T];
  real reward[N, T];

  int<lower=-1,upper=6> option1_t[N, T_t];
  int<lower=-1,upper=6> option2_t[N, T_t];
  int<lower=-1,upper=1> choice_t[N, T_t];
}

transformed data {
  // Default values to initialize the vector of expected values
  vector[6] initial_values;
  initial_values = rep_vector(0, 6);
}

parameters {
  // Group-level parameters
  vector[2] mu_pr;
  vector<lower=0>[2] sigma;

  // Subject-level parameters for Matt trick
  vector[N] alpha_pr;
  vector[N] beta_pr;
}

transformed parameters {
  vector<lower=0,upper=1>[N] alpha;
  vector<lower=0,upper=10>[N] beta;

  alpha = Phi_approx(mu_pr[1] + sigma[1] * alpha_pr);
  beta  = Phi_approx(mu_pr[2] + sigma[2] * beta_pr) * 10;
}

model {
  // Priors for group-level parameters
  mu_pr ~ normal(0, 1);
  sigma ~ normal(0, 0.2);

  // Priors for subject-level parameters
  alpha_pr ~ normal(0, 1);
  beta_pr  ~ normal(0, 1);

  for (i in 1:N) {
    int co;         // Chosen option
    real delta;     // Difference between two options (training)
    real delta_t;   // Difference between two options (test)
    real pe;        // Prediction error
    vector[6] ev;   // Expected values

    ev = initial_values;

    // Acquisition Phase
    for (t in 1:Tsubj[i]) {
      co = (choice[i, t] > 0) ? option1[i, t] : option2[i, t];

      // Luce choice rule
      delta = ev[option1[i, t]] - ev[option2[i, t]];
      choice[i, t] ~ bernoulli_logit(beta[i] * delta);

      pe = reward[i, t] - ev[co];
      ev[co] += alpha[i] * pe;
    }

    // Test phase
    for (u in 1:Tsubj_t[i]) {
      // Softmax
      delta_t = ev[option1_t[i, u]] - ev[option2_t[i, u]];
      choice_t[i, u] ~ bernoulli_logit(beta[i] * delta_t);
    }
  }
}

generated quantities {
  // For group-level parameters
  real<lower=0,upper=1>  mu_alpha;
  real<lower=0,upper=10> mu_beta;

  // For log-likelihood calculation
  real log_lik[N];

  // For posterior predictive check
  real y_pred[N, T_t];

 // Initialize all the variables to avoid NULL values
  for (i in 1:N) {
    for (t in 1:T) {
      y_pred[i, t]    = -1;
    }
  }

  mu_alpha = Phi_approx(mu_pr[1]);
  mu_beta  = Phi_approx(mu_pr[2]) * 10;

  {
    for (i in 1:N) {
      int co;         // Chosen option
      real delta;     // Difference between two options
      real delta_t;   // Difference between two options
      real pe;        // Prediction error
      vector[6] ev;   // Expected values

      ev = initial_values;
      log_lik[i] = 0;

      // Acquisition Phase
      for (t in 1:Tsubj[i]) {
        co = (choice[i, t] > 0) ? option1[i, t] : option2[i, t];

        // Luce choice rule
        delta = ev[option1[i, t]] - ev[option2[i, t]];
        log_lik[i] += bernoulli_logit_lpmf(choice[i, t] | beta[i] * delta);

        pe = reward[i, t] - ev[co];
        ev[co] += alpha[i] * pe;
      }

      // Test phase
      for (u in 1:Tsubj_t[i]) {
        delta_t = ev[option1_t[i, u]] - ev[option2_t[i, u]];
        // increment log-likelihood
        log_lik[i] += bernoulli_logit_lpmf(choice_t[i, u] | beta[i] * delta_t);
        // generate posterior prediction for current trial
        y_pred[i, u] =  bernoulli_logit_rng(beta[i] * delta_t);
      }
    }
  }
}
