// GLM for a Gaussian, Gamma, inverse Gaussian, or Beta outcome
functions {
#include /include/rsa_funcs.stan
  /**
   * Increments the log-posterior with the logarithm of a multivariate normal 
   * likelihood with a scalar standard deviation for all errors
   * Equivalent to normal_lpdf(y | intercept + X * beta + Z * b, sigma) but faster
   * @param coeff vector of coefficients (including intercept)
   * @param OLS precomputed vector of OLS coefficients (including intercept)
   * @param XtX precomputed matrix equal to crossprod(X) (including intercept)
   * @param SSR positive precomputed value of the sum of squared OLS residuals
   * @param sigma positive scalar for the standard deviation of the errors
   * @param N integer equal to the number of observations
   */
  real ll_mvn_ols(vector coeff, vector OLS, matrix XtX,
                  real SSR, real sigma, int N) {
    return -0.5 * (quad_form(XtX, coeff - OLS) + SSR) / square(sigma)
            - N * (log(sigma) + log(sqrt(2 * pi())));
  }

  /** 
  * test function for csr_matrix_times_vector
  *
  * @param m Integer number of rows
  * @param n Integer number of columns
  * @param w Vector (see reference manual)
  * @param v Integer array (see reference manual)
  * @param u Integer array (see reference manual)
  * @param b Vector that is multiplied from the left by the CSR matrix
  * @return A vector that is the product of the CSR matrix and b
  */
  vector test_csr_matrix_times_vector(int m, int n, vector w, 
                                      int[] v, int[] u, vector b) {
    return csr_matrix_times_vector(m, n, w, v, u, b); 
  }
}
data {
  // declares N, K, X, xbar, dense_X, nnz_x, w_x, v_x, u_x
  // dimensions
  int<lower=0> N;  // number of observations
  int<lower=0> K;  // number of predictors
  
  // data
  vector[K] xbar;               // predictor means
  int<lower=0,upper=1> dense_X; // flag for dense vs. sparse
  matrix[N,K] X[dense_X];       // centered predictor matrix in the dense case

  // stuff for the sparse case
  int<lower=0> nnz_X;                      // number of non-zero elements in the implicit X matrix
  vector[nnz_X] w_X;                       // non-zero elements in the implicit X matrix
  int<lower=0, upper = K - 1> v_X[nnz_X];  // column indices for w_X
  // where the non-zeros start in each row of X
  int<lower=0, upper = rows(w_X) + 1> u_X[dense_X ? 0 : N + 1]; 
  
  // smooths
  int<lower=0> K_smooth;
  matrix[N,K_smooth] S;
  int<lower=1> smooth_map[K_smooth];  int<lower=0> len_y;      // length of y
  real lb_y; // lower bound on y
  real<lower=lb_y> ub_y; // upper bound on y
  vector<lower=lb_y, upper=ub_y>[len_y] y; // continuous outcome
  int<lower=1,upper=4> family; // 1 gaussian, 2 gamma, 3 inv-gaussian, 4 beta
  // declares prior_PD, has_intercept, link, prior_dist, prior_dist_for_intercept
  // flag indicating whether to draw from the prior
  int<lower=0,upper=1> prior_PD;  // 1 = yes
  int<lower=0,upper=1> compute_mean_PPD; // 1 = yes
  
  // intercept
  int<lower=0,upper=1> has_intercept;  // 1 = yes
  
  // link function from location to linear predictor 
  int<lower=1> link;  // interpretation varies by .stan file
  
  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus, 
  //   5 = laplace, 6 = lasso, 7 = product_normal
  int<lower=0,upper=7> prior_dist;
  int<lower=0,upper=2> prior_dist_for_intercept;
  
  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = exponential
  int<lower=0,upper=3> prior_dist_for_aux;
  
  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = exponential
  int<lower=0,upper=3> prior_dist_for_smooth;
  // declares has_weights, weights, has_offset, offset
  // weights
  int<lower=0,upper=1> has_weights;  // 0 = No, 1 = Yes
  vector[has_weights ? N : 0] weights;
  
  // offset
  int<lower=0,upper=1> has_offset;  // 0 = No, 1 = Yes
  vector[has_offset ? N : 0] offset_;
  // declares prior_{mean, scale, df}, prior_{mean, scale, df}_for_intercept, prior_{mean, scale, df}_for_aux
  // hyperparameter values are set to 0 if there is no prior
  vector<lower=0>[K] prior_scale;
  real<lower=0> prior_scale_for_intercept;
  real<lower=0> prior_scale_for_aux;
  vector<lower=0>[K_smooth > 0 ? max(smooth_map) : 0] prior_scale_for_smooth;
  vector[K] prior_mean;
  real prior_mean_for_intercept;
  real<lower=0> prior_mean_for_aux;
  vector<lower=0>[K_smooth > 0 ? max(smooth_map) : 0] prior_mean_for_smooth;
  vector<lower=0>[K] prior_df;
  real<lower=0> prior_df_for_intercept;
  real<lower=0> prior_df_for_aux;
  vector<lower=0>[K_smooth > 0 ? max(smooth_map) : 0] prior_df_for_smooth;
  real<lower=0> global_prior_df;     // for hs priors only
  real<lower=0> global_prior_scale;  // for hs priors only
  real<lower=0> slab_df;     // for hs prior only
  real<lower=0> slab_scale;  // for hs prior only
  int<lower=2> num_normals[prior_dist == 7 ? K : 0];
  // declares t, p[t], l[t], q, len_theta_L, shape, scale, {len_}concentration, {len_}regularization
  // glmer stuff, see table 3 of
  // https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  int<lower=0> t;               // num. terms (maybe 0) with a | in the glmer formula
  int<lower=1> p[t];            // num. variables on the LHS of each |
  int<lower=1> l[t];            // num. levels for the factor(s) on the RHS of each |
  int<lower=0> q;               // conceptually equals \sum_{i=1}^t p_i \times l_i
  int<lower=0> len_theta_L;     // length of the theta_L vector

  // hyperparameters for glmer stuff; if t > 0 priors are mandatory
  vector<lower=0>[t] shape; 
  vector<lower=0>[t] scale;
  int<lower=0> len_concentration;
  real<lower=0> concentration[len_concentration];
  int<lower=0> len_regularization;
  real<lower=0> regularization[len_regularization];
  // declares num_not_zero, w, v, u
  int<lower=0> num_non_zero;  // number of non-zero elements in the Z matrix
  vector[num_non_zero] w;     // non-zero elements in the implicit Z matrix
  int<lower=0, upper=q-1> v[num_non_zero];               // column indices for w
  int<lower=0, upper=rows(w) + 1> u[t > 0 ? N + 1 : 0];  // where the non-zeros start in each row
  int<lower=0,upper=1> special_case;                     // is the only term (1|group)
  // betareg data
  int<lower=0, upper=1> has_intercept_z;  // presence of z intercept
  int<lower=0> link_phi;                  // link transformation for eta_z (0 => no z in model)
  int<lower=0> z_dim;                     // dimensions of z vars
  matrix[N, z_dim] betareg_z;             // matrix of z vars
  row_vector[z_dim] zbar;                 // mean of predictors
  // betareg hyperparameters
  int<lower=0,upper=7> prior_dist_z;
  int<lower=0,upper=2> prior_dist_for_intercept_z;
  vector<lower=0>[z_dim] prior_scale_z;
  real<lower=0> prior_scale_for_intercept_z;
  vector[z_dim] prior_mean_z;
  real prior_mean_for_intercept_z;
  vector<lower=0>[z_dim] prior_df_z;
  real<lower=0> prior_df_for_intercept_z;
  real<lower=0> global_prior_scale_z;
  real<lower=0> global_prior_df_z;
  real<lower=0> slab_df_z;
  real<lower=0> slab_scale_z;
  int<lower=2> num_normals_z[prior_dist_z == 7 ? z_dim : 0];
  int<lower=0,upper=10> SSfun; // nonlinear function indicator, 0 for identity
  vector[SSfun > 0  ? len_y : 0] input;
  vector[SSfun == 5 ? len_y : 0] Dose;
}
transformed data {
  vector[family == 3 ? len_y : 0] sqrt_y;
  vector[family == 3 ? len_y : 0] log_y;
  real sum_log_y = family == 1 ? not_a_number() : sum(log(y));
  int<lower=1> V[special_case ? t : 0, len_y] = make_V(len_y, special_case ? t : 0, v);
  int<lower=0> hs_z;                  // for tdata_betareg.stan
  int can_do_OLS = family == 1 && link == 1 && SSfun == 0 && has_offset == 0 && t == 0 && 
                   prior_PD == 0 && dense_X && N > 2 && len_y >= (has_intercept + K + K_smooth);
  vector[can_do_OLS ? has_intercept + K + K_smooth : 0] OLS;
  matrix[can_do_OLS ? has_intercept + K + K_smooth : 0, can_do_OLS ? has_intercept + K + K_smooth : 0] XtX;
  int can_do_normalidglm = K != 0 &&  // remove K!=0 after rstan includes this Stan bugfix: https://github.com/stan-dev/math/issues/1398 
                           can_do_OLS == 0 && family == 1 && link == 1 && 
                           SSfun == 0 && has_offset == 0 && dense_X && prior_PD == 0 && 
                           t == 0 && len_y < (has_intercept + K + K_smooth);
  matrix[can_do_normalidglm ? N : 0, can_do_normalidglm ? K + K_smooth : 0] XS;
  real SSR = not_a_number();
  // defines hs, len_z_T, len_var_group, delta, is_continuous, pos
  int<lower=0> len_z_T = 0;
  int<lower=0> len_var_group = sum(p) * (t > 0);
  int<lower=0> len_rho = sum(p) - t;
  int<lower=0, upper=1> is_continuous = 0; // changed in continuous.stan
  int<lower=1> pos = 1;
  real<lower=0> delta[len_concentration];
  int<lower=0> hs;
  if (prior_dist <= 2) hs = 0;
  else if (prior_dist == 3) hs = 2;
  else if (prior_dist == 4) hs = 4;
  else hs = 0;
  
  for (i in 1:t) {
    if (p[i] > 1) {
      for (j in 1:p[i]) {
        delta[pos] = concentration[j];
        pos += 1;
      }
    }
    for (j in 3:p[i]) len_z_T += p[i] - 1;
  }
  // defines hs_z
  if (prior_dist_z <= 2) hs_z = 0;
  else if (prior_dist_z == 3) hs_z = 2;
  else if (prior_dist_z == 4) hs_z = 4;
  else hs_z = 0;
  is_continuous = 1;

  if (family == 3) {
    sqrt_y = sqrt(y);
    log_y = log(y);
  }
  if (can_do_OLS) {
    matrix[N, has_intercept + K + K_smooth ] X_ = has_intercept ? append_col(rep_vector(1.0, N), 
                                                  (K_smooth > 0 ? append_col(X[1], S) : X[1])) : 
                                                  (K_smooth > 0 ? append_col(X[1], S) : X[1]);
    matrix[cols(X_), cols(X_)] R = qr_thin_R(X_);
    if (tail(diagonal(R), 1)[1] > 1e-16) {
      matrix[N, cols(R)] Q = qr_thin_Q(X_);
      XtX = crossprod(X_);
      OLS = mdivide_right_tri_low(y' * Q, R')';
      SSR = dot_self(y - X_ * OLS);
    } else can_do_OLS = 0;
  }
  if (can_do_normalidglm) {
    XS = K_smooth > 0 ? append_col(X[1], S) : X[1];
  }
}
parameters {
  real<lower=make_lower(family, link),upper=make_upper(family,link)> gamma[has_intercept];
  // declares z_beta, global, local, z_b, z_T, rho, zeta, tau
  vector[prior_dist == 7 ? sum(num_normals) : K] z_beta;
  vector[K_smooth] z_beta_smooth;
  vector<lower=0>[K_smooth > 0 ? smooth_map[K_smooth] : 0] smooth_sd_raw;
  real<lower=0> global[hs];
  vector<lower=0>[K] local[hs];
  real<lower=0> caux[hs > 0];
  vector<lower=0>[K] mix[prior_dist == 5 || prior_dist == 6];
  real<lower=0> one_over_lambda[prior_dist == 6];
  vector[q] z_b;
  vector[len_z_T] z_T;
  vector<lower=0,upper=1>[len_rho] rho;
  vector<lower=0>[len_concentration] zeta;
  vector<lower=0>[t] tau;
  real<lower=0> aux_unscaled; // interpretation depends on family!
  vector[prior_dist_z == 7 ? sum(num_normals_z) : z_dim] z_omega; // betareg z variable coefficients
  real<lower=(link_phi <= 1 ? negative_infinity() : 0)> gamma_z[has_intercept_z];  // betareg intercept
  real<lower=0> global_z[hs_z];
  vector<lower=0>[z_dim] local_z[hs_z];
  real<lower=0> caux_z[hs_z > 0];
  vector<lower=0>[z_dim] S_z[prior_dist_z == 5 || prior_dist_z == 6];
  real<lower=0> one_over_lambda_z[prior_dist_z == 6];
}
transformed parameters {
  // aux has to be defined first in the hs case
  real aux = prior_dist_for_aux == 0 ? aux_unscaled : (prior_dist_for_aux <= 2 ? 
             prior_scale_for_aux * aux_unscaled + prior_mean_for_aux :
             prior_scale_for_aux * aux_unscaled);

  vector[z_dim] omega; // used in tparameters_betareg.stan
  // defines beta, b, theta_L
  vector[K] beta;
  vector[K_smooth] beta_smooth;
  vector[K_smooth > 0 ? smooth_map[K_smooth] : 0] smooth_sd;
  vector[q] b;
  vector[len_theta_L] theta_L;
  if      (prior_dist == 0) beta = z_beta;
  else if (prior_dist == 1) beta = z_beta .* prior_scale + prior_mean;
  else if (prior_dist == 2) for (k in 1:K) {
    beta[k] = CFt(z_beta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
  }
  else if (prior_dist == 3) {
    real c2 = square(slab_scale) * caux[1];
    if (is_continuous == 1 && family == 1)
      beta = hs_prior(z_beta, global, local, global_prior_scale, aux, c2);
    else beta = hs_prior(z_beta, global, local, global_prior_scale, 1, c2);
  }
  else if (prior_dist == 4) {
    real c2 = square(slab_scale) * caux[1];
    if (is_continuous == 1 && family == 1)
      beta = hsplus_prior(z_beta, global, local, global_prior_scale, aux, c2);
    else beta = hsplus_prior(z_beta, global, local, global_prior_scale, 1, c2);
  }
  else if (prior_dist == 5) // laplace
    beta = prior_mean + prior_scale .* sqrt(2 * mix[1]) .* z_beta;
  else if (prior_dist == 6) // lasso
    beta = prior_mean + one_over_lambda[1] * prior_scale .* sqrt(2 * mix[1]) .* z_beta;
  else if (prior_dist == 7) { // product_normal
    int z_pos = 1;
    for (k in 1:K) {
      beta[k] = z_beta[z_pos];
      z_pos += 1;
      for (n in 2:num_normals[k]) {
        beta[k] *= z_beta[z_pos];
        z_pos += 1;
      }
      beta[k] *= prior_scale[k] ^ num_normals[k];
      beta[k] += prior_mean[k];
    }
  }

  if (K_smooth) {
    smooth_sd = prior_mean_for_smooth + prior_scale_for_smooth .* smooth_sd_raw;
    if (is_continuous && family == 1) smooth_sd *= aux;
    beta_smooth = z_beta_smooth .* smooth_sd[smooth_map];
  }
  if (prior_dist_z == 0) omega = z_omega;
  else if (prior_dist_z == 1) omega = z_omega .* prior_scale_z + prior_mean_z;
  else if (prior_dist_z == 2) for (k in 1:z_dim) {
    real left = CFt(omega[k], prior_df_z[k]); 
    omega[k] = left * prior_scale_z[k] + prior_mean_z[k];
  }
  else if (prior_dist_z == 3) 
    omega = hs_prior(z_omega, global_z, local_z, global_prior_scale, 
                     1, square(slab_scale_z) * caux_z[1]);
  else if (prior_dist_z == 4) 
    omega = hsplus_prior(z_omega, global_z, local_z, global_prior_scale, 1,
                         square(slab_scale_z) * caux_z[1]);
  else if (prior_dist_z == 5)
    omega = prior_mean_z + prior_scale_z .* sqrt(2 * S_z[1]) .* z_omega;
  else if (prior_dist_z == 6)
    omega = prior_mean_z + one_over_lambda_z[1] * prior_scale_z .* sqrt(2 * S_z[1]) .* z_omega;
  else if (prior_dist_z == 7) {
    int z_pos = 1;
    for (k in 1:z_dim) {
      omega[k] = z_omega[z_pos];
      z_pos += 1;
      for (n in 2:num_normals_z[k]) {
        omega[k] *= z_omega[z_pos];
        z_pos += 1;
      }
      omega[k] *= prior_scale_z[k] ^ num_normals_z[k];
      omega[k] += prior_mean_z[k];
    }
  }  
  if (prior_dist_for_aux == 0) // none
    aux = aux_unscaled;
  else {
    aux = prior_scale_for_aux * aux_unscaled;
    if (prior_dist_for_aux <= 2) // normal or student_t
      aux += prior_mean_for_aux;
  }

  if (t > 0) {
    if (special_case == 1) {
      int start = 1;
      theta_L = scale .* tau * aux;
      if (t == 1) b = theta_L[1] * z_b;
      else for (i in 1:t) {
        int end = start + l[i] - 1;
        b[start:end] = theta_L[i] * z_b[start:end];
        start = end + 1;
      }
    }
    else {
      theta_L = make_theta_L(len_theta_L, p, 
                             aux, tau, scale, zeta, rho, z_T);
      b = make_b(z_b, theta_L, p, l);
    }
  }
}
model {
  if (can_do_OLS) {
    vector[cols(XtX)] coeff = has_intercept ? append_row(to_vector(gamma), 
                                              (K_smooth > 0 ? append_row(beta, beta_smooth) : beta)) : 
                                              (K_smooth > 0 ? append_row(beta, beta_smooth) : beta);
    target += ll_mvn_ols(coeff, OLS, XtX, SSR, aux, N);
  } else if (can_do_normalidglm) {
    vector[K + K_smooth] coeff = K_smooth > 0 ? append_row(beta, beta_smooth) : beta;
    target += normal_id_glm_lpdf(y | XS, has_intercept ? gamma[1] : 0.0, coeff, aux);
  } else if (prior_PD == 0) {
    vector[link_phi > 0 ? N : 0] eta_z; // beta regression linear predictor for phi
  vector[N] eta;  // linear predictor
  if (K > 0) {
    if (dense_X) eta = X[1] * beta;
    else eta = csr_matrix_times_vector(N, K, w_X, v_X, u_X, beta);
  }
  else eta = rep_vector(0.0, N);
  if (has_offset == 1) eta += offset_;
  if (K_smooth) eta += S * beta_smooth;
    if (t > 0) {
			if (special_case) for (i in 1:t) eta += b[V[i]];
			else eta += csr_matrix_times_vector(N, q, w, v, u, b);
    }
    if (has_intercept == 1) {
      if ((family == 1 || link == 2) || (family == 4 && link != 5)) eta += gamma[1];
      else if (family == 4 && link == 5) eta += gamma[1] - max(eta);
      else eta += gamma[1] - min(eta);
    }
    else {
		// correction to eta if model has no intercept (because X is centered)
			eta += dot_product(xbar, beta); 
    }
    if (SSfun > 0) { // nlmer
      matrix[len_y, K] P = reshape_vec(eta, len_y, K);
      if (SSfun < 5) {
        if (SSfun <= 2) {
          if (SSfun == 1) target += normal_lpdf(y | SS_asymp(input, P), aux);
          else target += normal_lpdf(y | SS_asympOff(input, P), aux);
        }
        else if (SSfun == 3) target += normal_lpdf(y | SS_asympOrig(input, P), aux);
        else {
          for (i in 1:len_y) P[i,1] += exp(P[i,3]); // ordering constraint
          target += normal_lpdf(y | SS_biexp(input, P), aux);
        }
      }
      else {
        if (SSfun <= 7) {
          if (SSfun == 5) target += normal_lpdf(y | SS_fol(Dose, input, P), aux);
          else if (SSfun == 6) target += normal_lpdf(y | SS_fpl(input, P), aux);
          else target += normal_lpdf(y | SS_gompertz(input, P), aux);
        }
        else {
          if (SSfun == 8) target += normal_lpdf(y | SS_logis(input, P), aux);
          else if (SSfun == 9) target += normal_lpdf(y | SS_micmen(input, P), aux);
          else target += normal_lpdf(y | SS_weibull(input, P), aux);
        }
      }
    }
    else if (has_weights == 0) { // unweighted log-likelihoods
			if (family == 4 && z_dim > 0 && link_phi > 0) {
				eta_z = betareg_z * omega;
			}
			else if (family == 4 && z_dim == 0 && has_intercept_z == 1){
				eta_z = rep_vector(0.0, N); 
			}
      // adjust eta_z according to links
      if (has_intercept_z == 1) {
        if (link_phi > 1) {
          eta_z += gamma_z[1] - min(eta_z);
        }
        else {
          eta_z += gamma_z[1];
        }
      }
      else { // has_intercept_z == 0
				if (link_phi > 1) {
					eta_z += dot_product(zbar, omega) - min(eta_z);
				}
				else {
					eta_z += dot_product(zbar, omega);
				}
      }
      if (family == 1) {
        if (link == 1) 
          target += normal_lpdf(y | eta, aux);
        else if (link == 2) 
          target += normal_lpdf(y | exp(eta), aux);
        else 
          target += normal_lpdf(y | inv(eta), aux);
      }
      else if (family == 2) {
        target += GammaReg(y, eta, aux, link, sum_log_y);
      }
      else if (family == 3) {
        target += inv_gaussian(y, linkinv_inv_gaussian(eta, link), 
                               aux, sum_log_y, sqrt_y);
      }
      else if (family == 4 && link_phi == 0) {
        vector[N] mu;
        mu = linkinv_beta(eta, link);
        target += beta_lpdf(y | mu * aux, (1 - mu) * aux);
      }
      else if (family == 4 && link_phi > 0) {
        vector[N] mu;
        vector[N] mu_z;
        mu = linkinv_beta(eta, link);
        mu_z = linkinv_beta_z(eta_z, link_phi);
        target += beta_lpdf(y | rows_dot_product(mu, mu_z), 
                            rows_dot_product((1 - mu) , mu_z));
      }
    }
    else { // weighted log-likelihoods
      vector[N] summands;
      if (family == 1) summands = pw_gauss(y, eta, aux, link);
      else if (family == 2) summands = pw_gamma(y, eta, aux, link);
      else if (family == 3) summands = pw_inv_gaussian(y, eta, aux, link, log_y, sqrt_y);
      else if (family == 4 && link_phi == 0) summands = pw_beta(y, eta, aux, link);
      else if (family == 4 && link_phi > 0) summands = pw_beta_z(y, eta, eta_z, link, link_phi);
      target += dot_product(weights, summands);
    }
  }

  // Log-priors
  if (prior_dist_for_aux > 0 && prior_scale_for_aux > 0) {
    real log_half = -0.693147180559945286;
    if (prior_dist_for_aux == 1)
      target += normal_lpdf(aux_unscaled | 0, 1) - log_half;
    else if (prior_dist_for_aux == 2)
      target += student_t_lpdf(aux_unscaled | prior_df_for_aux, 0, 1) - log_half;
    else 
     target += exponential_lpdf(aux_unscaled | 1);
  }
    
  // Log-priors for coefficients
       if (prior_dist == 1) target += normal_lpdf(z_beta | 0, 1);
  else if (prior_dist == 2) target += normal_lpdf(z_beta | 0, 1); // Student t via Cornish-Fisher expansion
  else if (prior_dist == 3) { // hs
    real log_half = -0.693147180559945286;
    target += normal_lpdf(z_beta | 0, 1);
    target += normal_lpdf(local[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
    target += normal_lpdf(global[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
    target += inv_gamma_lpdf(caux | 0.5 * slab_df, 0.5 * slab_df);
  }
  else if (prior_dist == 4) { // hs+
    real log_half = -0.693147180559945286;
    target += normal_lpdf(z_beta | 0, 1);
    target += normal_lpdf(local[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
    target += normal_lpdf(local[3] | 0, 1) - log_half;
    // unorthodox useage of prior_scale as another df hyperparameter
    target += inv_gamma_lpdf(local[4] | 0.5 * prior_scale, 0.5 * prior_scale);
    target += normal_lpdf(global[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
    target += inv_gamma_lpdf(caux | 0.5 * slab_df, 0.5 * slab_df);
  }
  else if (prior_dist == 5) { // laplace
    target += normal_lpdf(z_beta | 0, 1);
    target += exponential_lpdf(mix[1] | 1);
  }
  else if (prior_dist == 6) { // lasso
    target += normal_lpdf(z_beta | 0, 1);
    target += exponential_lpdf(mix[1] | 1);
    target += chi_square_lpdf(one_over_lambda[1] | prior_df[1]);
  }
  else if (prior_dist == 7) { // product_normal
    target += normal_lpdf(z_beta | 0, 1);
  }
  /* else prior_dist is 0 and nothing is added */
  
  // Log-prior for intercept  
  if (has_intercept == 1) {
    if (prior_dist_for_intercept == 1)  // normal
      target += normal_lpdf(gamma | prior_mean_for_intercept, prior_scale_for_intercept);
    else if (prior_dist_for_intercept == 2)  // student_t
      target += student_t_lpdf(gamma | prior_df_for_intercept, prior_mean_for_intercept, 
                               prior_scale_for_intercept);
    /* else prior_dist is 0 and nothing is added */
  }

  if (K_smooth) {
    target += normal_lpdf(z_beta_smooth | 0, 1);
    if (prior_dist_for_smooth > 0) {
      real log_half = -0.693147180559945286;
      if (prior_dist_for_smooth == 1) 
        target += normal_lpdf(smooth_sd_raw | 0, 1) - log_half;
      else if (prior_dist_for_smooth == 2)
        target += student_t_lpdf(smooth_sd_raw | prior_df_for_smooth, 0, 1) - log_half;
      else if (prior_dist_for_smooth == 3) 
        target += exponential_lpdf(smooth_sd_raw | 1);
    }
  }
  // Log-priors for coefficients
  if (prior_dist_z == 1)  target += normal_lpdf(z_omega | 0, 1);
  else if (prior_dist_z == 2) target += normal_lpdf(z_omega | 0, 1);
  else if (prior_dist_z == 3) { // hs
    real log_half = -0.693147180559945286;
    target += normal_lpdf(z_omega | 0, 1);
    target += normal_lpdf(local_z[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(local_z[2] | 0.5 * prior_df_z, 0.5 * prior_df_z);
    target += normal_lpdf(global_z[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(global_z[2] | 0.5 * global_prior_df_z, 0.5 * global_prior_df_z);
    target += inv_gamma_lpdf(caux_z | 0.5 * slab_df_z, 0.5 * slab_df_z);
  }
  else if (prior_dist_z == 4) { // hs+
    real log_half = -0.693147180559945286;
    target += normal_lpdf(z_omega | 0, 1);
    target += normal_lpdf(local_z[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(local_z[2] | 0.5 * prior_df_z, 0.5 * prior_df_z);
    target += normal_lpdf(local_z[3] | 0, 1) - log_half;
    // unorthodox useage of prior_scale as another df hyperparameter
    target += inv_gamma_lpdf(local_z[4] | 0.5 * prior_scale_z, 0.5 * prior_scale_z);
    target += normal_lpdf(global_z[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(global_z[2] | 0.5, 0.5);
    target += inv_gamma_lpdf(caux_z | 0.5 * slab_df_z, 0.5 * slab_df_z);
  }
  else if (prior_dist_z == 5) { // laplace
    target += normal_lpdf(z_omega | 0, 1);
    target += exponential_lpdf(S_z[1] | 1);
  }
  else if (prior_dist_z == 6) { // lasso
    target += normal_lpdf(z_omega | 0, 1);
    target += exponential_lpdf(S_z[1] | 1);
    target += chi_square_lpdf(one_over_lambda_z[1] | prior_df_z[1]);
  }
  else if (prior_dist_z == 7) { // product_normal
    target += normal_lpdf(z_omega | 0, 1);
  }
  /* else prior_dist is 0 and nothing is added */
  
  // Log-prior for intercept  
  if (has_intercept_z == 1) {
    if (prior_dist_for_intercept_z == 1)  // normal
      target += normal_lpdf(gamma_z | prior_mean_for_intercept_z, prior_scale_for_intercept_z);
    else if (prior_dist_for_intercept_z == 2)  // student_t
      target += student_t_lpdf(gamma_z | prior_df_for_intercept_z, prior_mean_for_intercept_z, 
                               prior_scale_for_intercept_z);
    /* else prior_dist is 0 and nothing is added */
  }
  if (t > 0) {
    real dummy = decov_lp(z_b, z_T, rho, zeta, tau, 
                          regularization, delta, shape, t, p);
  }
}
generated quantities {
  real mean_PPD = compute_mean_PPD ? 0 : negative_infinity();
  real alpha[has_intercept];
  real omega_int[has_intercept_z];
  
  if (has_intercept == 1) {
    if (dense_X) alpha[1] = gamma[1] - dot_product(xbar, beta);
    else alpha[1] = gamma[1];
  }
  if (has_intercept_z == 1) { 
    omega_int[1] = gamma_z[1] - dot_product(zbar, omega);  // adjust betareg intercept 
  }
  
  if (compute_mean_PPD) {
    vector[N] eta_z;
		vector[N] eta;  // linear predictor
		if (K > 0) {
			if (dense_X) eta = X[1] * beta;
			else eta = csr_matrix_times_vector(N, K, w_X, v_X, u_X, beta);
		}
		else eta = rep_vector(0.0, N);
		if (has_offset == 1) eta += offset_;
		if (K_smooth) eta += S * beta_smooth;
    if (t > 0) {
			if (special_case) for (i in 1:t) eta += b[V[i]];
			else eta += csr_matrix_times_vector(N, q, w, v, u, b);
    }
    if (has_intercept == 1) {
      if (make_lower(family,link) == negative_infinity() &&
          make_upper(family,link) == positive_infinity()) eta += gamma[1];
      else if (family == 4 && link == 5) {
        real max_eta = max(eta);
        alpha[1] -= max_eta;
        eta += gamma[1] - max_eta;
      }
      else {
        real min_eta = min(eta);
        alpha[1] -= min_eta;
        eta += gamma[1] - min_eta;
      }
    }
    else {
			if (link_phi > 1) {
				eta_z += dot_product(zbar, omega) - min(eta_z);
			}
			else {
				eta_z += dot_product(zbar, omega);
			}
    }
	  if (family == 4 && z_dim > 0 && link_phi > 0) {
			eta_z = betareg_z * omega;
		}
		else if (family == 4 && z_dim == 0 && has_intercept_z == 1){
			eta_z = rep_vector(0.0, N); 
		}    
		// adjust eta_z according to links
    if (has_intercept_z == 1) {
      if (link_phi > 1) {
        omega_int[1] -= min(eta_z);
        eta_z += gamma_z[1] - min(eta_z);
      }
      else {
        eta_z += gamma_z[1];
      }
    }
    else { // has_intercept_z == 0
			if (link_phi > 1) {
				eta_z += dot_product(zbar, omega) - min(eta_z);
			}
			else {
				eta_z += dot_product(zbar, omega);
			}    
		}    
    if (SSfun > 0) { // nlmer
      vector[len_y] eta_nlmer;
      matrix[len_y, K] P;      
      P = reshape_vec(eta, len_y, K);
      if (SSfun < 5) {
        if (SSfun <= 2) {
          if (SSfun == 1) eta_nlmer = SS_asymp(input, P);
          else eta_nlmer = SS_asympOff(input, P);
        }
        else if (SSfun == 3) eta_nlmer = SS_asympOrig(input, P);
        else eta_nlmer = SS_biexp(input, P);
      }
      else {
        if (SSfun <= 7) {
          if (SSfun == 5) eta_nlmer = SS_fol(Dose, input, P);
          else if (SSfun == 6) eta_nlmer = SS_fpl(input, P);
          else eta_nlmer = SS_gompertz(input, P);
        }
        else {
          if (SSfun == 8) eta_nlmer = SS_logis(input, P);
          else if (SSfun == 9) eta_nlmer = SS_micmen(input, P);
          else eta_nlmer = SS_weibull(input, P);
        }
      }
      for (n in 1:len_y) mean_PPD += normal_rng(eta_nlmer[n], aux);
    }
    else if (family == 1) {
      vector[N] mu = link > 1 ? linkinv_gauss(eta, link) : eta;
      for (n in 1:len_y) mean_PPD += normal_rng(mu[n], aux);
    }
    else if (family == 2) {
      vector[N] mu = link > 1 ? linkinv_gamma(eta, link) : eta;
      for (n in 1:len_y) mean_PPD += gamma_rng(aux, aux / mu[n]);
    }
    else if (family == 3) {
      vector[N] mu = link > 1 ? linkinv_inv_gaussian(eta, link) : eta;
      for (n in 1:len_y) mean_PPD += inv_gaussian_rng(mu[n], aux);
    }
    else if (family == 4 && link_phi == 0) { 
      vector[N] mu = linkinv_beta(eta, link);
      for (n in 1:N) {
        real mu_n = mu[n];
        if (aux <= 0) mean_PPD += bernoulli_rng(0.5);
        else if (mu_n >= 1) mean_PPD += 1;
        else if (mu_n > 0)
          mean_PPD += beta_rng(mu_n * aux, (1 - mu_n) * aux);
      }
    }
    else if (family == 4 && link_phi > 0) {
      vector[N] mu = linkinv_beta(eta, link);
      vector[N] phi = linkinv_beta_z(eta_z, link_phi);
      for (n in 1:N) {
        real mu_n = mu[n];
        real aux_n = phi[n];
        if (aux_n <= 0) mean_PPD += bernoulli_rng(0.5);
        else if (mu_n >= 1) mean_PPD += 1;
        else if (mu_n > 0)
          mean_PPD += beta_rng(mu_n * aux_n, (1 - mu_n) * aux_n);
      }
    }
    mean_PPD /= len_y;
  }
}
