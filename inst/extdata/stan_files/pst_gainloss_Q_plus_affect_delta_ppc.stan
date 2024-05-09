// -----------------------------------------------------------------------------
// Gain-loss Q-learning model for PST training data + change in affect rating
//------------------------------------------------------------------------------
// References:
//// https://www.pnas.org/doi/10.1073/pnas.1407535111
//// https://osf.io/9g2zw (replication of the above with Stan code)
//// https://www.nature.com/articles/s41562-023-01519-7 (mood drift over time)
//------------------------------------------------------------------------------

data {
  int<lower=1> N, T, I;               // # participants, max # of trials, max (theoretical) # of intervening trials
  array[N] int Tsubj;                 // # of trials for acquisition phase

  array[N, T] int option1;            // LHS option (1-6)
  array[N, T] int option2;            // RHS option (1-6)
  array[N, T] int choice;             // choice (1 = chose option 1, 3, or 5)
  matrix[N, T] reward;                // coded as 1 (reward) or -1 (no reward)

  array[N] row_vector[T] affect;      // includes 0 and 1, needs to be transformed
  array[N, T] int question;           // from 1 to 3 (happy, confident, engaged)
  array[N] row_vector[T] ovl_time;    // in hours to keep weights relatively small
  array[N, T] int int_trials;         // intervening trials between ratings
}

transformed data {
  // default values to initialize the vectors/matrices of EVs/PEs
  vector[6] inits;
  row_vector[T] neg_ones;
  inits = rep_vector(0, 6);
  neg_ones = rep_row_vector(-1, T);

  // transform affect to be strictly between 0 and 1 (Smith & Verkuilen, 2006)
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
  sigma_ql ~ normal(0, 0.2);

  // hyperpriors on the weights
  for (q in 1:3) {
    mu_wt[q]     ~ normal(0, 1);
    sigma_wt[q]  ~ exponential(0.1);
  }

  // hyperpriors on the beta distribution precision
  aff_mu_phi    ~ normal(0, 1);
  aff_sigma_phi ~ exponential(0.1);

  // priors on QL parameters
  alpha_pos_pr ~ normal(0, 1);
  alpha_neg_pr ~ normal(0, 1);
  beta_pr      ~ normal(0, 1);

  // priors on the mean-level weights and beta distribution precision
  for (q in 1:3) {
    w0_pr[:, q]       ~ normal(0, 1);
    w1_o_pr[:, q]     ~ normal(0, 1);
    w2_i_pr[:, q]     ~ normal(0, 1);
    w3_i_pr[:, q]     ~ normal(0, 1);
    sigma_w2_pr[:, q] ~ exponential(0.1);
    sigma_w3_pr[:, q] ~ exponential(0.1);
    phi_pr[:, q]      ~ normal(0, 1);
    for (j in 1:I) {
      w2_pr[:, q, j]  ~ normal(0, 1);
      w3_pr[:, q, j]  ~ normal(0, 1);
    }
  }

  for (i in 1:N) {
    int ti;
    ti = Tsubj[i];           // Number of trials for participant

    int co;                  // Chosen option
    int qn;                  // Question number
    real pe;                 // Prediction error
    real alpha;              // Learning rate (positive or negative)
    
    real eta;                // Linear predictor for affect at trial t
    vector[6] ev;            // Expected values per symbol
    vector[ti] delta;        // Difference in EVs between options

    vector[ti] ev_vec;       // Vector of summed EVs by trial
    vector[ti] pe_vec;       // Vector of summed PEs by trial
    vector[ti] phi_vec;      // Vector of beta distribution precision by trial
    vector[ti] eta_vec;      // Linear predictor for affect by trial

    vector[ti] aff_mu_cond;  // Conditional mean of the beta distribution
    vector[ti] shape_a;      // Beta distribution shape parameter alpha
    vector[ti] shape_b;      // Beta distribution shape parameter beta

    ev = inits;

    // initialise at machine precision to ensure shape parameters > 0
    shape_a = rep_vector(machine_precision(), ti);
    shape_b = rep_vector(machine_precision(), ti);

    // calculate relevant quantities for each trial
    for (t in 1:ti) {
      co = (choice[i, t] > 0) ? option1[i, t] : option2[i, t];
      qn = question[i, t];

      // Luce choice rule (EVs of unseen options assumed not to affect choice)
      delta[t] = ev[option1[i, t]] - ev[option2[i, t]];
      
      pe = reward[i, t] - ev[co];    
      alpha = (pe >= 0) ? alpha_pos[i] : alpha_neg[i];
      ev[co] += alpha * pe;

      // store summed EVs and PEs for this trial
      ev_vec[t] = ev[co];
      pe_vec[t] = pe;

      // initial value of linear predictor for affect at trial t
      eta = w0[i, question[i, t]] + w1_o[i, question[i, t]] * ovl_time[i, t];

      // add time-dependent coefficients to linear predictor depending on the number of intervening trials
      for (j in 1:int_trials[i, t]) {
        int t1 = t + 1; // so as to include the current trial
        eta += (
          w2[i][question[i, t1 - j], j] * ev_vec[t1 - j] + 
          w3[i][question[i, t1 - j], j] * pe_vec[t1 - j]
        );
      }

      // store AR scale parameters and affect estimates for this trial
      eta_vec[t] = eta;
      phi_vec[t] = phi[i, question[i, t]];
    }

    // increment log density for choice for participant i
    choice[i, :ti] ~ bernoulli_logit(beta[i] * delta);

    // calculate conditional mean of the beta distribution
    aff_mu_cond = inv_logit(eta_vec);

    // calculate beta distribution shape parameters (Smith & Verkuilen, 2006)
    shape_a += aff_mu_cond .* phi_vec;
    shape_b += phi_vec .* (1 - aff_mu_cond);
    
    // increment log density for affect for participant i
    affect_tr[i][:ti] ~ beta(shape_a, shape_b);
  }
}

generated quantities {
  // group-level parameter means
  real<lower=0,upper=1>  mu_alpha_pos;
  real<lower=0,upper=1>  mu_alpha_neg;
  real<lower=0,upper=10> mu_beta;

  vector[3] mu_w0;
  vector[3] mu_w1_o;

  // initialise log-likelihood vector and posterior prediction matrix
  vector[N] log_lik;
  array[N] row_vector[T] y_pred;

  y_pred = rep_array(neg_ones, N);

  // calculate moments of the group-level posterior distributions
  mu_alpha_pos = Phi_approx(mu_ql[1]);
  mu_alpha_neg = Phi_approx(mu_ql[2]);
  mu_beta      = Phi_approx(mu_ql[3]) * 10;

  mu_w0    = mu_wt[:, 1];
  mu_w1_o  = mu_wt[:, 2];

  {
    for (i in 1:N) {
      int ti;
      ti = Tsubj[i];           // Number of trials for participant

      int co;                  // Chosen option
      int qn;                  // Question number
      real pe;                 // Prediction error
      real alpha;              // Learning rate (positive or negative)
      
      real eta;                // Linear predictor for affect at trial t
      vector[6] ev;            // Expected values per symbol
      vector[ti] delta;        // Difference in EVs between options

      vector[ti] ev_vec;       // Vector of summed EVs by trial
      vector[ti] pe_vec;       // Vector of summed PEs by trial
      vector[ti] phi_vec;      // Vector of beta distribution precision by trial
      vector[ti] eta_vec;      // Linear predictor for affect by trial

      vector[ti] aff_mu_cond;  // Conditional mean of the beta distribution
      vector[ti] shape_a;      // Beta distribution shape parameter alpha
      vector[ti] shape_b;      // Beta distribution shape parameter beta

      ev = inits;

      // initialise at machine precision to ensure shape parameters > 0
      shape_a = rep_vector(machine_precision(), ti);
      shape_b = rep_vector(machine_precision(), ti);

      log_lik[i] = 0;

      // calculate relevant quantities for each trial
      for (t in 1:ti) {
        co = (choice[i, t] > 0) ? option1[i, t] : option2[i, t];
        qn = question[i, t];

        // Luce choice rule (EVs of unseen options assumed not to affect choice)
        delta[t] = ev[option1[i, t]] - ev[option2[i, t]];
        
        pe = reward[i, t] - ev[co];    
        alpha = (pe >= 0) ? alpha_pos[i] : alpha_neg[i];
        ev[co] += alpha * pe;

        // store summed EVs and PEs for this trial
        ev_vec[t] = ev[co];
        pe_vec[t] = pe;

        // initial value of linear predictor for affect at trial t
        eta = w0[i, question[i, t]] + w1_o[i, question[i, t]] * ovl_time[i, t];

        // add time-dependent coefficients to linear predictor
        for (j in 1:int_trials[i, t]) {
          int t1 = t + 1; // so as to include the current trial
          eta += (
            w2[i][question[i, t1 - j], j] * ev_vec[t1 - j] + 
            w3[i][question[i, t1 - j], j] * pe_vec[t1 - j]
          );
        }

        // store AR scale parameters and affect estimates for this trial
        eta_vec[t] = eta;
        phi_vec[t] = phi[i, question[i, t]];
      }

      // increment log likelihood for choice for participant i
      log_lik[i] += bernoulli_logit_lpmf(choice[i, :ti] | beta[i] * delta);

      // calculate conditional mean of the beta distribution
      aff_mu_cond = inv_logit(eta_vec);

      // calculate beta distribution shape parameters (Smith & Verkuilen, 2006)
      shape_a += aff_mu_cond .* phi_vec;
      shape_b += phi_vec .* (1 - aff_mu_cond);
      
      // increment log likelihood for affect for participant i
      log_lik[i] += beta_lpdf(affect_tr[i][:ti] | shape_a, shape_b);

      // generate posterior predictions on original scale
      y_pred[i][:ti] =
        ((to_row_vector(beta_rng(shape_a, shape_b)) * N) - 0.5) / (N - 1);
    }
  }
}
