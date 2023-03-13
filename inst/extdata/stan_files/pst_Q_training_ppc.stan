// Q-learning model for PST training data (incl. posterior predictive checks)
data {
  int<lower=1> N, T;          // # participants, max # of trials
  array[N] int Tsubj;         // # of trials for acquisition phase

  array[N, T] int option1;    // LHS option (1-6)
  array[N, T] int option2;    // RHS option (1-6)
  array[N, T] int choice;     // choice (1 = chose option 1, 3, or 5)
  matrix[N, T] reward;        // coded as 1 (reward) or -1 (no reward)
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
    real delta;     // Difference between two options
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
  }
}

generated quantities {
  // For group-level parameters
  real<lower=0,upper=1>  mu_alpha;
  real<lower=0,upper=10> mu_beta;

  // initialise log-likelihood vector and posterior prediction matrix
  vector[N] log_lik;
  vector[N] neg_ones;
  matrix[N, T] y_pred;

  neg_ones = rep_vector(-1, N);
  y_pred   = rep_matrix(neg_ones, T);

  mu_alpha = Phi_approx(mu_pr[1]);
  mu_beta  = Phi_approx(mu_pr[2]) * 10;

  {
    for (i in 1:N) {
      int co;         // Chosen option
      real delta;     // Difference between two options
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

        // generate posterior prediction for current trial
        y_pred[i, t] =  bernoulli_logit_rng(beta[i] * delta);

        pe = reward[i, t] - ev[co];
        ev[co] += alpha[i] * pe;
      }
    }
  }
}
