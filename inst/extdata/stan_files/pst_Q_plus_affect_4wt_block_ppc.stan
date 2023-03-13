// -----------------------------------------------------------------------------
// Q-learning model for PST training data + affect + block number
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
  array[N, T] int block_no;   // from 1 to 6
}

transformed data {
  // default values to initialize the vectors/matrices of EVs/PEs
  vector[6] inits;
  row_vector[T] zeros;
  inits = rep_vector(0, 6);
  zeros = rep_row_vector(0, T);

  int s;
  matrix[N, T] affect_tr;
  s = N*T;

  // transform affect to be strictly between 0 and 1 (Smith & Verkuilen, 2006)
  for (t in 1:T) {
    affect_tr[:, t] = ((affect[:, t]*(s-1))+0.5)/s;
  }
}

parameters {
  // group-level RL parameters
  vector[2] mu_ql;
  vector<lower=0>[2] sigma_ql;

  // group-level weights
  matrix[4, 3] mu_wt;
  matrix<lower=0>[4, 3] sigma_wt;

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
  matrix[N, 3] w2;
  matrix[N, 3] w3;
  matrix<lower=0, upper=1>[N, 3] gamma;
  matrix[N, 3] phi;

  for (p in 1:3) {
    w0[:, p]    = mu_wt[1, p] + sigma_wt[1, p] * w0_pr[:, p];
    w1_o[:, p]  = mu_wt[2, p] + sigma_wt[2, p] * w1_o_pr[:, p];
    w2[:, p]    = mu_wt[3, p] + sigma_wt[3, p] * w2_pr[:, p];
    w3[:, p]    = mu_wt[4, p] + sigma_wt[4, p] * w3_pr[:, p];
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
    w2_pr[:, p]   ~ normal(0, 1);
    w3_pr[:, p]   ~ normal(0, 1);
    gm_pr[:, p]   ~ normal(0, 1);
    phi_pr[:, p]  ~ normal(0, 1);
  }

  for (i in 1:N) {
    int ti;
    ti = Tsubj[i];           // Number of trials for participant i

    int co;                  // Chosen option
    real pe;                 // Prediction error

    vector[6] ev;            // Expected values per symbol
    vector[ti] delta;        // Difference in EVs between options

    vector[ti] w0_vec;       // Vector of w0 weights
    vector[ti] w1_o_vec;     // Vector of w1_o weights
    vector[ti] w2_vec;       // Vector of w2 weights
    vector[ti] w3_vec;       // Vector of w3 weights
    vector[ti] phi_vec;      // Vector of beta distribution precision
    
    row_vector[ti] z_vec;    // Vector of zeros
    row_vector[ti] ev_vec;   // Vector of expected values
    row_vector[ti] pe_vec;   // Vector of prediction errors
    
    matrix[ti, 3] dcy_mat;   // Weighting of previous trials, by question

    vector[ti] ev_dcy;       // Vector of summed decayed EVs by trial
    vector[ti] pe_dcy;       // Vector of summed decayed PEs by trial

    vector[ti] aff_mu_cond;  // Conditional mean of the beta distribution
    vector[ti] shape_a;      // Beta distribution shape parameter alpha
    vector[ti] shape_b;      // Beta distribution shape parameter beta

    z_vec   = zeros[:ti];
    ev      = inits;
    ev_vec  = z_vec;
    pe_vec  = z_vec;
    dcy_mat = rep_matrix(z_vec', 3);
    ev_dcy  = z_vec';
    pe_dcy  = z_vec';

    // initialise at machine precision to ensure shape parameters > 0
    shape_a = rep_vector(machine_precision(), ti);
    shape_b = rep_vector(machine_precision(), ti);

    // calculate relevant quantities for each trial
    for (t in 1:ti) {
      co = (choice[i, t] > 0) ? option1[i, t] : option2[i, t];

      // Luce choice rule (i.e., EVs of non-seen options assumed not to matter)
      delta[t] = ev[option1[i, t]] - ev[option2[i, t]];

      pe = reward[i, t] - ev[co];
      pe_vec[t] = pe;

      ev[co] += alpha[i] * pe;
      ev_vec[t] = ev[co];

      // store decay factor weighting for each question (gammas are different)
      for (q in 1:3) {
        dcy_mat[t, q] = pow(gamma[i, q], t - 1);
      }
      
      // store weights and beta distribution precision for convenience
      w0_vec[t]  = w0[i, question[i, t]];
      w1_o_vec[t] = w1_o[i, question[i, t]];
      w2_vec[t]  = w2[i, question[i, t]];
      w3_vec[t]  = w3[i, question[i, t]];
      phi_vec[t] = phi[i, question[i, t]];

      // store decayed EVs and PEs (i.e., gamma weighted sum over prev. trials)
      ev_dcy[t] = reverse(ev_vec[:t]) * dcy_mat[:t, question[i, t]];
      pe_dcy[t] = reverse(pe_vec[:t]) * dcy_mat[:t, question[i, t]];
    }

    // increment log density for choice for participant i
    choice[i, :ti] ~ bernoulli_logit(beta[i] * delta);
    
    // calculate conditional mean of the beta distribution
    aff_mu_cond = inv_logit(
      w0_vec + 
      w1_o_vec .* to_vector(block_no[i]) + 
      (w2_vec .* ev_dcy) + 
      (w3_vec .* pe_dcy)
    );

    // calculate beta distribution shape parameters (Smith & Verkuilen, 2006)
    shape_a += aff_mu_cond .* phi_vec;
    shape_b += phi_vec .* (1 - aff_mu_cond);
    
    // increment log density for affect for participant i
    affect_tr[i, :ti] ~ beta(shape_a, shape_b);
  }
}

generated quantities {
  // group-level parameter means
  real<lower=0,upper=1>  mu_alpha;
  real<lower=0,upper=10> mu_beta;

  vector[3] mu_w0;
  vector[3] mu_w1_o;
  vector[3] mu_w2;
  vector[3] mu_w3;
  vector[3] mu_gamma;

  // initialise log-likelihood vector and posterior prediction matrix
  vector[N] log_lik;
  vector[N] neg_ones;
  matrix[N, T] y_pred;

  neg_ones = rep_vector(-1, N);
  y_pred = rep_matrix(neg_ones, T);

  // calculate moments of the group-level posterior distributions
  mu_alpha = Phi_approx(mu_ql[1]);
  mu_beta  = Phi_approx(mu_ql[2]) * 10;

  for (p in 1:3) {
    mu_w0[p]    = mu_wt[1, p];
    mu_w1_o[p]  = mu_wt[2, p];
    mu_w2[p]    = mu_wt[3, p];
    mu_w3[p]    = mu_wt[4, p];
    mu_gamma[p] = Phi_approx(mu_gm[p]);
  }

  // calculate log-likelihoods and posterior predictions
  for (i in 1:N) {
    int ti;
    ti = Tsubj[i];           // Number of trials for participant i

    int co;                  // Chosen option
    real pe;                 // Prediction error

    vector[6] ev;            // Expected values per symbol
    vector[ti] delta;        // Difference in EVs between options

    vector[ti] w0_vec;       // Vector of w0 weights
    vector[ti] w1_o_vec;     // Vector of w1_o weights
    vector[ti] w2_vec;       // Vector of w2 weights
    vector[ti] w3_vec;       // Vector of w3 weights
    vector[ti] phi_vec;      // Vector of beta distribution precision
    
    row_vector[ti] z_vec;    // Vector of zeros
    row_vector[ti] ev_vec;   // Vector of expected values
    row_vector[ti] pe_vec;   // Vector of prediction errors
    
    matrix[ti, 3] dcy_mat;   // Weighting of previous trials, by question

    vector[ti] ev_dcy;       // Vector of summed decayed EVs by trial
    vector[ti] pe_dcy;       // Vector of summed decayed PEs by trial

    vector[ti] aff_mu_cond;  // Conditional mean of the beta distribution
    vector[ti] shape_a;      // Beta distribution shape parameter alpha
    vector[ti] shape_b;      // Beta distribution shape parameter beta

    z_vec   = zeros[:ti];
    ev      = inits;
    ev_vec  = z_vec;
    pe_vec  = z_vec;
    dcy_mat = rep_matrix(z_vec', 3);
    ev_dcy  = z_vec';
    pe_dcy  = z_vec';
    
    // initialise at machine precision to ensure shape parameters > 0
    shape_a = rep_vector(machine_precision(), ti);
    shape_b = rep_vector(machine_precision(), ti);

    log_lik[i] = 0;

    // calculate relevant quantities for each trial
    for (t in 1:ti) {
      co = (choice[i, t] > 0) ? option1[i, t] : option2[i, t];

      // Luce choice rule (i.e., EVs of non-seen options assumed not to matter)
      delta[t] = ev[option1[i, t]] - ev[option2[i, t]];
      
      pe = reward[i, t] - ev[co];
      pe_vec[t] = pe;

      ev[co] += alpha[i] * pe;
      ev_vec[t] = ev[co];

      // store decay factor weighting for each question (gammas are different)
      for (q in 1:3) {
        dcy_mat[t, q] = pow(gamma[i, q], t - 1);
      }
      
      // store weights and beta distribution precision for convenience
      w0_vec[t]  = w0[i, question[i, t]];
      w1_o_vec[t] = w1_o[i, question[i, t]];
      w2_vec[t]  = w2[i, question[i, t]];
      w3_vec[t]  = w3[i, question[i, t]];
      phi_vec[t] = phi[i, question[i, t]];

      // store decayed EVs and PEs (i.e., gamma weighted sum over prev. trials)
      ev_dcy[t] = reverse(ev_vec[:t]) * dcy_mat[:t, question[i, t]];
      pe_dcy[t] = reverse(pe_vec[:t]) * dcy_mat[:t, question[i, t]];
    }
    
    // increment log likelihood for choice for participant i
    log_lik[i] += bernoulli_logit_lpmf(choice[i, :ti] | beta[i] * delta);

    // calculate conditional mean of the beta distribution
    aff_mu_cond = inv_logit(
      w0_vec + 
      w1_o_vec .* to_vector(block_no[i]) + 
      (w2_vec .* ev_dcy) + 
      (w3_vec .* pe_dcy)
    );

    // calculate beta distribution shape parameters (Smith & Verkuilen, 2006)
    shape_a += aff_mu_cond .* phi_vec;
    shape_b += phi_vec .* (1 - aff_mu_cond);
    
    // increment log likelihood for affect for participant i
    log_lik[i] += beta_lpdf(affect_tr[i, :ti] | shape_a, shape_b);

    // generate posterior predictions
    y_pred[i, :ti] = to_row_vector(beta_rng(shape_a, shape_b));
  }

  // backtransform predictions to be on the original scale for completeness
  for (t in 1:T) {
    y_pred[:, t] = ((y_pred[:, t]*s)-0.5)/(s-1);
  }
}
