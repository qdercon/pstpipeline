// Q-learning model for PST training data (for separate posterior prediction)
data {
  int<lower=1> N;             // Number of subjects
  int<lower=1> T;             // Maximum # of trials
  int<lower=1> Tsubj[N];      // # of trials for acquisition phase

  int<lower=-1,upper=6> option1[N, T];
  int<lower=-1,upper=6> option2[N, T];
  int<lower=-1,upper=1> choice[N, T];
  real reward[N, T];
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

generated quantities {
  // For group-level parameters
  real<lower=0,upper=1>  mu_alpha;
  real<lower=0,upper=10> mu_beta;

  // For posterior predictive check
  real y_pred[N, T];

  // Initialize all the variables to avoid NULL values
  for (i in 1:N) {
    for (t in 1:T) {
      y_pred[i, t]    = -1;
    }
  }

  {
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
		
        // generate posterior prediction for current trial
        y_pred[i, t] =  bernoulli_logit_rng(beta[i] * delta);

        pe = reward[i, t] - ev[co];
        ev[co] += alpha[i] * pe;
      }
    }
  }
}