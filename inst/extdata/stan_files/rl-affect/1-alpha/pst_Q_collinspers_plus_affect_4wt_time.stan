// -----------------------------------------------------------------------------
// Q-learning model for PST training data + affect + overall time
// Collins-style perseveration: alpha_neg = alpha * (1 - pers)
//------------------------------------------------------------------------------
// References:
//// https://www.pnas.org/doi/10.1073/pnas.1407535111
//// https://osf.io/9g2zw (replication of the above with Stan code)
//// https://psyarxiv.com/v86bx (Beta distribution response)
//// https://psyarxiv.com/bwv58 (passage-of-time-dysphoria, inspired time pars)
//------------------------------------------------------------------------------

data {
  int<lower=0,upper=1> run_gq;      // 1 to run generated quantities, 0 to skip
  int<lower=1> N, T;                // # participants, max # of trials
  array[N] int Tsubj;               // # of trials for acquisition phase

  array[N, T] int option1;          // LHS option (1-6)
  array[N, T] int option2;          // RHS option (1-6)
  array[N, T] int choice;           // choice (1 = chose option 1, 3, or 5)
  matrix[N, T] reward;              // coded as 1 (reward) or -1 (no reward)

  array[N] row_vector[T] affect;    // includes 0 and 1, needs to be transformed
  array[N, T] int question;         // from 1 to 3 (happy, confident, engaged)
  array[N] row_vector[T] ovl_time;  // in hours to keep weights relatively small
}

transformed data {
  vector[6] inits;
  vector[T] zeros;
  row_vector[3] init_sum;
  row_vector[T] neg_ones;

  inits = rep_vector(0, 6);
  zeros = rep_vector(0, T);
  init_sum = rep_row_vector(0, 3);
  neg_ones = rep_row_vector(-1, T);

  array[N] row_vector[T] affect_tr;
  for (i in 1:N) {
    affect_tr[i] = ((affect[i] * (N - 1)) + 0.5) / N;
  }

  array[N] vector[T] ovl_time_tr;
  for (i in 1:N) {
    ovl_time_tr[i] = to_vector(ovl_time[i]);
  }
}

parameters {
  vector[3] mu_ql;
  vector<lower=0>[3] sigma_ql;

  matrix[3, 4] mu_wt;
  matrix<lower=0>[3, 4] sigma_wt;

  vector[3] mu_gm;
  vector<lower=0>[3] sigma_gm;

  vector[3] aff_mu_phi;
  vector<lower=0>[3] aff_sigma_phi;

  vector[N] alpha_pr;
  vector[N] pers_pr;
  vector[N] beta_pr;

  matrix[N, 3] w0_pr;
  matrix[N, 3] w1_o_pr;
  matrix[N, 3] w2_pr;
  matrix[N, 3] w3_pr;
  matrix[N, 3] gm_pr;

  matrix[N, 3] phi_pr;
}

transformed parameters {
  vector<lower=0,upper=1>[N] alpha;
  vector<lower=0,upper=1>[N] pers;
  vector<lower=0,upper=1>[N] alpha_neg;
  vector<lower=0,upper=10>[N] beta;

  alpha     = Phi_approx(mu_ql[1] + sigma_ql[1] * alpha_pr);
  pers      = Phi_approx(mu_ql[2] + sigma_ql[2] * pers_pr);
  alpha_neg = alpha .* (1 - pers);
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
  mu_ql    ~ normal(0, 1);
  sigma_ql ~ normal(0, 0.2);

  for (q in 1:3) {
    mu_wt[q]    ~ normal(0, 1);
    sigma_wt[q] ~ exponential(0.1);
  }

  mu_gm    ~ normal(0, 1);
  sigma_gm ~ exponential(0.1);

  aff_mu_phi    ~ normal(0, 1);
  aff_sigma_phi ~ exponential(0.1);

  alpha_pr ~ normal(0, 1);
  pers_pr  ~ normal(0, 1);
  beta_pr  ~ normal(0, 1);

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
    ti = Tsubj[i];

    int co;
    int qn;
    real pe;
    real alpha_t;

    vector[6] ev;
    vector[ti] delta;

    vector[ti] w0_vec;
    vector[ti] w1_o_vec;
    vector[ti] w2_vec;
    vector[ti] w3_vec;
    vector[ti] phi_vec;
    vector[ti] z_vec;

    row_vector[3] ev_sum;
    row_vector[3] pe_sum;

    vector[ti] ev_dcy;
    vector[ti] pe_dcy;

    vector[ti] aff_mu_cond;
    vector[ti] shape_a;
    vector[ti] shape_b;

    z_vec   = zeros[:ti];
    ev      = inits;
    ev_sum  = init_sum;
    pe_sum  = init_sum;
    ev_dcy  = z_vec;
    pe_dcy  = z_vec;

    shape_a = rep_vector(machine_precision(), ti);
    shape_b = rep_vector(machine_precision(), ti);

    for (t in 1:ti) {
      co = (choice[i, t] > 0) ? option1[i, t] : option2[i, t];
      qn = question[i, t];

      delta[t] = ev[option1[i, t]] - ev[option2[i, t]];

      pe = reward[i, t] - ev[co];
      pe_sum += pe;

      alpha_t = (pe >= 0) ? alpha[i] : alpha_neg[i];
      ev[co] += alpha_t * pe;
      ev_sum += ev[co];

      ev_dcy[t] = ev_sum[qn];
      pe_dcy[t] = pe_sum[qn];

      w0_vec[t]   = w0[i, qn];
      w1_o_vec[t] = w1_o[i, qn];
      w2_vec[t]   = w2[i, qn];
      w3_vec[t]   = w3[i, qn];
      phi_vec[t]  = phi[i, qn];

      ev_sum = ev_sum .* gamma[i, :];
      pe_sum = pe_sum .* gamma[i, :];
    }

    choice[i, :ti] ~ bernoulli_logit(beta[i] * delta);

    aff_mu_cond = inv_logit(
      w0_vec +
      w1_o_vec .* ovl_time_tr[i][:ti] +
      w2_vec .* ev_dcy +
      w3_vec .* pe_dcy
    );

    shape_a += aff_mu_cond .* phi_vec;
    shape_b += phi_vec .* (1 - aff_mu_cond);

    affect_tr[i][:ti] ~ beta(shape_a, shape_b);
  }
}

generated quantities {
  real<lower=0,upper=1>  mu_alpha;
  real<lower=0,upper=1>  mu_pers;
  real<lower=0,upper=1>  mu_alpha_neg;
  real<lower=0,upper=10> mu_beta;

  vector[3] mu_w0;
  vector[3] mu_w1_o;
  vector[3] mu_w2;
  vector[3] mu_w3;
  vector[3] mu_gamma;

  vector[N] log_lik;
  array[N] row_vector[T] y_pred;

  y_pred = rep_array(neg_ones, N);

  mu_alpha     = Phi_approx(mu_ql[1]);
  mu_pers      = Phi_approx(mu_ql[2]);
  mu_alpha_neg = mu_alpha * (1 - mu_pers);
  mu_beta      = Phi_approx(mu_ql[3]) * 10;

  mu_w0    = mu_wt[:, 1];
  mu_w1_o  = mu_wt[:, 2];
  mu_w2    = mu_wt[:, 3];
  mu_w3    = mu_wt[:, 4];
  mu_gamma = Phi_approx(mu_gm);

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

  for (i in 1:N) {
    int ti;
    ti = Tsubj[i];

    int co;
    int qn;
    real pe;
    real alpha_t;

    vector[6] ev;
    vector[ti] delta;

    vector[ti] w0_vec;
    vector[ti] w1_o_vec;
    vector[ti] w2_vec;
    vector[ti] w3_vec;
    vector[ti] phi_vec;
    vector[ti] z_vec;

    row_vector[3] ev_sum;
    row_vector[3] pe_sum;

    vector[ti] ev_dcy;
    vector[ti] pe_dcy;

    vector[ti] aff_mu_cond;
    vector[ti] shape_a;
    vector[ti] shape_b;

    z_vec   = zeros[:ti];
    ev      = inits;
    ev_sum  = init_sum;
    pe_sum  = init_sum;
    ev_dcy  = z_vec;
    pe_dcy  = z_vec;

    shape_a = rep_vector(machine_precision(), ti);
    shape_b = rep_vector(machine_precision(), ti);

    log_lik[i] = 0;

    for (t in 1:ti) {
      co = (choice[i, t] > 0) ? option1[i, t] : option2[i, t];
      qn = question[i, t];

      delta[t] = ev[option1[i, t]] - ev[option2[i, t]];

      pe = reward[i, t] - ev[co];
      pe_sum += pe;

      alpha_t = (pe >= 0) ? alpha[i] : alpha_neg[i];
      ev[co] += alpha_t * pe;
      ev_sum += ev[co];

      ev_dcy[t] = ev_sum[qn];
      pe_dcy[t] = pe_sum[qn];

      w0_vec[t]   = w0[i, qn];
      w1_o_vec[t] = w1_o[i, qn];
      w2_vec[t]   = w2[i, qn];
      w3_vec[t]   = w3[i, qn];
      phi_vec[t]  = phi[i, qn];

      ev_sum = ev_sum .* gamma[i, :];
      pe_sum = pe_sum .* gamma[i, :];
    }

    log_lik[i] += bernoulli_logit_lpmf(choice[i, :ti] | beta[i] * delta);

    aff_mu_cond = inv_logit(
      w0_vec +
      w1_o_vec .* ovl_time_tr[i][:ti] +
      w2_vec .* ev_dcy +
      w3_vec .* pe_dcy
    );

    shape_a += aff_mu_cond .* phi_vec;
    shape_b += phi_vec .* (1 - aff_mu_cond);

    log_lik[i] += beta_lpdf(affect_tr[i][:ti] | shape_a, shape_b);

    y_pred[i][:ti] =
      ((to_row_vector(beta_rng(shape_a, shape_b)) * N) - 0.5) / (N - 1);
  }
}
