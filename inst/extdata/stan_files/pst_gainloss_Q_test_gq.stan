// Gain-loss Q-learning model for PST test data (for separate posterior prediction)
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
  vector[3] mu_pr;
  vector<lower=0>[3] sigma;

  // Subject-level parameters for Matt trick
  vector[N] alpha_pos_pr;
  vector[N] alpha_neg_pr;
  vector[N] beta_pr;
}

transformed parameters {
  vector<lower=0,upper=1>[N] alpha_pos;
  vector<lower=0,upper=1>[N] alpha_neg;
  vector<lower=0,upper=10>[N] beta;

  alpha_pos = Phi_approx(mu_pr[1] + sigma[1] * alpha_pos_pr);
  alpha_neg = Phi_approx(mu_pr[2] + sigma[2] * alpha_neg_pr);
  beta      = Phi_approx(mu_pr[3] + sigma[3] * beta_pr) * 10;
}

generated quantities {
  // initialise posterior prediction matrix
  vector[N] neg_ones;
  matrix[N, T_t] y_pred;

  neg_ones = rep_vector(-1, N);
  y_pred   = rep_matrix(neg_ones, T_t);
  
  {
    for (i in 1:N) {
      int co;         // Chosen option
      real delta;     // Difference between two options
      real delta_t;   // Difference between two options
      real pe;        // Prediction error
      real alpha;
      vector[6] ev;   // Expected values

      ev = initial_values;

      // Acquisition Phase
      for (t in 1:Tsubj[i]) {
        co = (choice[i, t] > 0) ? option1[i, t] : option2[i, t];

        // Luce choice rule
        delta = ev[option1[i, t]] - ev[option2[i, t]];

        pe = reward[i, t] - ev[co];
        alpha = (pe >= 0) ? alpha_pos[i] : alpha_neg[i];
        ev[co] += alpha * pe;
      }

      // Test phase
      for (u in 1:Tsubj_t[i]) {
        delta_t = ev[option1_t[i, u]] - ev[option2_t[i, u]];
		
        // generate posterior prediction for current trial
        y_pred[i, u] =  bernoulli_logit_rng(beta[i] * delta_t);
      }
    }
  }
}
