// Q-learning model for PST test data
data {
  int<lower=1> N, T, T_t;       // # participants, max # trials, max # test trials
  array[N] int Tsubj;           // # of trials for acquisition phase
  array[N] int Tsubj_t;         // # of trials for test phase

  array[N, T] int option1;      // LHS option (1-6)
  array[N, T] int option2;      // RHS option (1-6)
  array[N, T] int choice;       // choice (1 = chose option 1, 3, or 5)
  matrix[N, T] reward;          // coded as 1 (reward) or -1 (no reward)

  array[N, T_t] int option1_t;  // LHS option (1-6) in test phase
  array[N, T_t] int option2_t;  // RHS option (1-6) in test phase
  array[N, T_t] int choice_t;   // choice (1 = chose option 1, 3, or 5) in test
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

  // initialise log-likelihood vector
  vector[N] log_lik;

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
      }
    }
  }
}
