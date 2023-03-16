// -----------------------------------------------------------------------------
// Gain-loss Q-learning model for PST training data + affect + overall time
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
}

transformed data {
  // default values to initialize the vectors/matrices of EVs/PEs
  vector[6] inits;
  row_vector[T] zeros;
  inits = rep_vector(0, 6);
  zeros = rep_row_vector(0, T);

  // transform affect to be strictly between 0 and 1 (Smith & Verkuilen, 2006)
  matrix[N, T] affect_tr;
  affect_tr = ((affect * (N - 1)) + 0.5) / N;
}

parameters {
  // group-level RL parameters
  vector[3] mu_ql;
  vector<lower=0>[3] sigma_ql;

  // group-level weights
  matrix[3, 4] mu_wt; // 3 questions x 4 weights
  matrix<lower=0>[3, 4] sigma_wt;

  // group-level parameters for decay factor (gamma)
  vector[3] mu_gm;
  vector<lower=0>[3] sigma_gm;

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

  alpha_pos = Phi_approx(mu_ql[1] + sigma_ql[1] * alpha_pos_pr);
  alpha_neg = Phi_approx(mu_ql[2] + sigma_ql[2] * alpha_neg_pr);
  beta      = Phi_approx(mu_ql[3] + sigma_ql[3] * beta_pr) * 10;

  matrix[N, 3] w0;
  matrix[N, 3] w1_o;
  matrix[N, 3] w2;
  matrix[N, 3] w3;
  matrix<lower=0, upper=1>[N, 3] gamma;
  matrix[N, 3] phi;

  for (q in 1:3) {
    w0[:, q]    = mu_wt[q, 1] + sigma_wt[q, 1] * w0_pr[:, q];
    w1_o[:, q]  = mu_wt[q, 2] + sigma_wt[q, 2] * w1_o_pr[:, q];
    w2[:, q]    = mu_wt[q, 3] + sigma_wt[q, 3] * w2_pr[:, q];
    w3[:, q]    = mu_wt[q, 4] + sigma_wt[q, 4] * w3_pr[:, q];
    gamma[:, q] = Phi_approx(mu_gm[q] + sigma_gm[q] * gm_pr[:, q]);
    phi[:, q]   = exp(aff_mu_phi[q] + aff_sigma_phi[q] * phi_pr[:, q]);
  }
}

model {
  // hyperpriors on QL parameters
  mu_ql    ~ normal(0, 1);
  sigma_ql ~ normal(0, 0.2);

  // hyperpriors on the weights
  for (q in 1:3) {
    mu_wt[q]    ~ normal(0, 1);
    sigma_wt[q] ~ exponential(0.1);
  }

  // hyperpriors on gamma
  mu_gm    ~ normal(0, 1);
  sigma_gm ~ exponential(0.1);

  // hyperpriors on the beta distribution precision
  aff_mu_phi    ~ normal(0, 1);
  aff_sigma_phi ~ exponential(0.1);

  // priors on QL parameters
  alpha_pos_pr ~ normal(0, 1);
  alpha_neg_pr ~ normal(0, 1);
  beta_pr      ~ normal(0, 1);

  // priors on the weights + gamma + beta distribution precision
  for (q in 1:3) {
    w0_pr[:, q]   ~ normal(0, 1);
    w1_o_pr[:, q] ~ normal(0, 1);
    w2_pr[:, q]   ~ normal(0, 1);
    w3_pr[:, q]   ~ normal(0, 1);
    gm_pr[:, q]   ~ normal(0, 1);
    phi_pr[:, q]  ~ normal(0, 1);
  }

  for (i in 1:N) {
    int ti;
    ti = Tsubj[i];           // Number of trials for participant i

    int co;                  // Chosen option
    int qn;                  // Question number
    real pe;                 // Prediction error
    real alpha;              // Learning rate (positive or negative)
    
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
    
    matrix[ti, 3] dcy_qn;   // Weighting of previous trials, by question

    vector[ti] ev_dcy;       // Vector of summed decayed EVs by trial
    vector[ti] pe_dcy;       // Vector of summed decayed PEs by trial

    vector[ti] aff_mu_cond;  // Conditional mean of the beta distribution
    vector[ti] shape_a;      // Beta distribution shape parameter alpha
    vector[ti] shape_b;      // Beta distribution shape parameter beta

    z_vec   = zeros[:ti];
    ev      = inits;
    ev_vec  = z_vec;
    pe_vec  = z_vec;
    ev_dcy  = z_vec';
    pe_dcy  = z_vec';

    // compute decay matrix for each question
    for (q in 1:3) {
      for (t in 1:ti) {
        dcy_qn[t, q] = pow(gamma[i, q], t-1);
      }
    }
    
    // initialise at machine precision to ensure shape parameters > 0
    shape_a = rep_vector(machine_precision(), ti);
    shape_b = rep_vector(machine_precision(), ti);

    // calculate relevant quantities for each trial
    for (t in 1:ti) {
      co = (choice[i, t] > 0) ? option1[i, t] : option2[i, t];
      qn = question[i, t];

      // Luce choice rule (i.e., EVs of non-seen options assumed not to matter)
      delta[t] = ev[option1[i, t]] - ev[option2[i, t]];
      
      pe = reward[i, t] - ev[co];
      pe_vec[t] = pe;

      alpha = (pe >= 0) ? alpha_pos[i] : alpha_neg[i];
      ev[co] += alpha * pe;
      ev_vec[t] = ev[co];
      
      // store weights and beta distribution precision for convenience
      w0_vec[t]   = w0[i, qn];
      w1_o_vec[t] = w1_o[i, qn];
      w2_vec[t]   = w2[i, qn];
      w3_vec[t]   = w3[i, qn];
      phi_vec[t]  = phi[i, qn];

      // store decayed EVs and PEs (i.e., gamma weighted sum over prev. trials)
      ev_dcy[t] = dot_product(reverse(ev_vec[:t]), dcy_qn[:t, qn]);
      pe_dcy[t] = dot_product(reverse(pe_vec[:t]), dcy_qn[:t, qn]);
    }
    
    // increment log density for choice for participant i
    choice[i, :ti] ~ bernoulli_logit(beta[i] * delta);
    
    // calculate conditional mean of the beta distribution
    aff_mu_cond = inv_logit(
      w0_vec + 
      w1_o_vec .* to_vector(ovl_time[i, :ti]) + 
      w2_vec .* ev_dcy +
      w3_vec .* pe_dcy
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
  real<lower=0,upper=1>  mu_alpha_pos;
  real<lower=0,upper=1>  mu_alpha_neg;
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
  mu_alpha_pos = Phi_approx(mu_ql[1]);
  mu_alpha_neg = Phi_approx(mu_ql[2]);
  mu_beta      = Phi_approx(mu_ql[3]) * 10;

  mu_w0    = mu_wt[:, 1];
  mu_w1_o  = mu_wt[:, 2];
  mu_w2    = mu_wt[:, 3];
  mu_w3    = mu_wt[:, 4];
  mu_gamma = Phi_approx(mu_gm);

  // difference in weights between questions
  vector[3] w0_diff;
  vector[3] w1_o_diff;
  vector[3] w2_diff;
  vector[3] w3_diff;
  vector[3] gamma_diff;

  array[3] int bsl = { 1, 1, 2 };
  array[3] int comp = { 2, 3, 3 };

  w0_diff    = mu_w0[bsl] - mu_w0[comp];
  w1_o_diff  = mu_w1_o[bsl] - mu_w1_o[comp];
  w2_diff    = mu_w2[bsl] - mu_w2[comp];
  w3_diff    = mu_w3[bsl] - mu_w3[comp];
  gamma_diff = mu_gamma[bsl] - mu_gamma[comp];

  // calculate log-likelihoods and posterior predictions
  for (i in 1:N) {
    int ti;
    ti = Tsubj[i];           // Number of trials for participant i

    int co;                  // Chosen option
    int qn;                  // Question number
    real pe;                 // Prediction error
    real alpha;              // Learning rate (positive or negative)
    
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
    
    matrix[ti, 3] dcy_qn;   // Weighting of previous trials, by question

    vector[ti] ev_dcy;       // Vector of summed decayed EVs by trial
    vector[ti] pe_dcy;       // Vector of summed decayed PEs by trial

    vector[ti] aff_mu_cond;  // Conditional mean of the beta distribution
    vector[ti] shape_a;      // Beta distribution shape parameter alpha
    vector[ti] shape_b;      // Beta distribution shape parameter beta

    z_vec   = zeros[:ti];
    ev      = inits;
    ev_vec  = z_vec;
    pe_vec  = z_vec;
    ev_dcy  = z_vec';
    pe_dcy  = z_vec';

    // compute decay matrix for each question
    for (q in 1:3) {
      for (t in 1:ti) {
        dcy_qn[t, q] = pow(gamma[i, q], t-1);
      }
    }
    
    // initialise at machine precision to ensure shape parameters > 0
    shape_a = rep_vector(machine_precision(), ti);
    shape_b = rep_vector(machine_precision(), ti);

    log_lik[i] = 0;

    // calculate relevant quantities for each trial
    for (t in 1:ti) {
      co = (choice[i, t] > 0) ? option1[i, t] : option2[i, t];
      qn = question[i, t];

      // Luce choice rule (i.e., EVs of non-seen options assumed not to matter)
      delta[t] = ev[option1[i, t]] - ev[option2[i, t]];
      
      pe = reward[i, t] - ev[co];
      pe_vec[t] = pe;

      alpha = (pe >= 0) ? alpha_pos[i] : alpha_neg[i];
      ev[co] += alpha * pe;
      ev_vec[t] = ev[co];
      
      // store weights and beta distribution precision for convenience
      w0_vec[t]   = w0[i, qn];
      w1_o_vec[t] = w1_o[i, qn];
      w2_vec[t]   = w2[i, qn];
      w3_vec[t]   = w3[i, qn];
      phi_vec[t]  = phi[i, qn];

      // store decayed EVs and PEs (i.e., gamma weighted sum over prev. trials)
      ev_dcy[t] = dot_product(reverse(ev_vec[:t]), dcy_qn[:t, qn]);
      pe_dcy[t] = dot_product(reverse(pe_vec[:t]), dcy_qn[:t, qn]);
    }

    // increment log likelihood for choice for participant i
    log_lik[i] += bernoulli_logit_lpmf(choice[i, :ti] | beta[i] * delta);
    
    // calculate conditional mean of the beta distribution
    aff_mu_cond = inv_logit(
      w0_vec + 
      w1_o_vec .* to_vector(ovl_time[i, :ti]) + 
      w2_vec .* ev_dcy + 
      w3_vec .* pe_dcy
    );

    // calculate beta distribution shape parameters (Smith & Verkuilen, 2006)
    shape_a += aff_mu_cond .* phi_vec;
    shape_b += phi_vec .* (1 - aff_mu_cond);
    
    // increment log likelihood for affect for participant i
    log_lik[i] += beta_lpdf(affect_tr[i, :ti] | shape_a, shape_b);

    // generate posterior predictions
    y_pred[i, :ti] = to_row_vector(beta_rng(shape_a, shape_b));
  }

  // backtransform predictions to be on the original scale
  y_pred = ((y_pred * N) - 0.5) / (N - 1);
}
