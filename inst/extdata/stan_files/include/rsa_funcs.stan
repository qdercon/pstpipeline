	/* for multiple .stan files */
  
  /** 
   * Create group-specific block-diagonal Cholesky factor, see section 2 of
   * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
   * @param len_theta_L An integer indicating the length of returned vector, 
   *   which lme4 denotes as m
   * @param p An integer array with the number variables on the LHS of each |
   * @param dispersion Scalar standard deviation of the errors, calles sigma by lme4
   * @param tau Vector of scale parameters whose squares are proportional to the 
   *   traces of the relative covariance matrices of the group-specific terms
   * @param scale Vector of prior scales that are multiplied by elements of tau
   * @param zeta Vector of positive parameters that are normalized into simplexes
   *   and multiplied by the trace of the covariance matrix to produce variances
   * @param rho Vector of radii in the onion method for creating Cholesky factors
   * @param z_T Vector used in the onion method for creating Cholesky factors
   * @return A vector that corresponds to theta in lme4
   */
  vector make_theta_L(int len_theta_L, int[] p, real dispersion,
                      vector tau, vector scale, vector zeta,
                      vector rho, vector z_T) {
    vector[len_theta_L] theta_L;
    int zeta_mark = 1;
    int rho_mark = 1;
    int z_T_mark = 1;
    int theta_L_mark = 1;

    // each of these is a diagonal block of the implicit Cholesky factor
    for (i in 1:size(p)) { 
      int nc = p[i];
      if (nc == 1) { // "block" is just a standard deviation
        theta_L[theta_L_mark] = tau[i] * scale[i] * dispersion;
        // unlike lme4, theta[theta_L_mark] includes the dispersion term in it
        theta_L_mark += 1;
      }
      else { // block is lower-triangular               
        matrix[nc,nc] T_i; 
        real std_dev;
        real T21;
        real trace_T_i = square(tau[i] * scale[i] * dispersion) * nc;
        vector[nc] pi = segment(zeta, zeta_mark, nc); // gamma(zeta | shape, 1)
        pi /= sum(pi);                            // thus dirichlet(pi | shape)
        
        // unlike lme4, T_i includes the dispersion term in it
        zeta_mark += nc;
        std_dev = sqrt(pi[1] * trace_T_i);
        T_i[1,1] = std_dev;
        
        // Put a correlation into T_i[2,1] and scale by std_dev
        std_dev = sqrt(pi[2] * trace_T_i);
        T21 = 2.0 * rho[rho_mark] - 1.0;
        rho_mark += 1;
        T_i[2,2] = std_dev * sqrt(1.0 - square(T21));
        T_i[2,1] = std_dev * T21;
        
        for (r in 2:(nc - 1)) { // scaled onion method to fill T_i
          int rp1 = r + 1;
          vector[r] T_row = segment(z_T, z_T_mark, r);
          real scale_factor = sqrt(rho[rho_mark] / dot_self(T_row)) * std_dev;
          z_T_mark += r;
          std_dev = sqrt(pi[rp1] * trace_T_i);
          for(c in 1:r) T_i[rp1,c] = T_row[c] * scale_factor;
          T_i[rp1,rp1] = sqrt(1.0 - rho[rho_mark]) * std_dev;
          rho_mark += 1;
        }
        
        // now vech T_i
        for (c in 1:nc) for (r in c:nc) {
          theta_L[theta_L_mark] = T_i[r,c];
          theta_L_mark += 1;
        }
      }
    }
    return theta_L;
  }
  
  /** 
  * Create group-specific coefficients, see section 2 of
  * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  *
  * @param z_b Vector whose elements are iid normal(0,sigma) a priori
  * @param theta Vector with covariance parameters as defined in lme4
  * @param p An integer array with the number variables on the LHS of each |
  * @param l An integer array with the number of levels for the factor(s) on 
  *   the RHS of each |
  * @return A vector of group-specific coefficients
  */
  vector make_b(vector z_b, vector theta_L, int[] p, int[] l) {
    vector[rows(z_b)] b;
    int b_mark = 1;
    int theta_L_mark = 1;
    for (i in 1:size(p)) {
      int nc = p[i];
      if (nc == 1) {
        real theta_L_start = theta_L[theta_L_mark];
        for (s in b_mark:(b_mark + l[i] - 1)) 
          b[s] = theta_L_start * z_b[s];
        b_mark += l[i];
        theta_L_mark += 1;
      }
      else {
        matrix[nc,nc] T_i = rep_matrix(0, nc, nc);
        for (c in 1:nc) {
          T_i[c,c] = theta_L[theta_L_mark];
          theta_L_mark += 1;
          for(r in (c+1):nc) {
            T_i[r,c] = theta_L[theta_L_mark];
            theta_L_mark += 1;
          }
        }
        for (j in 1:l[i]) {
          vector[nc] temp = T_i * segment(z_b, b_mark, nc);
          b_mark -= 1;
          for (s in 1:nc) b[b_mark + s] = temp[s];
          b_mark += nc + 1;
        }
      }
    }
    return b;
  }

  /** 
   * Prior on group-specific parameters
   *
   * @param z_b A vector of primitive coefficients
   * @param z_T A vector of primitives for the unit vectors in the onion method
   * @param rho A vector radii for the onion method
   * @param zeta A vector of primitives for the simplexes
   * @param tau A vector of scale parameters
   * @param regularization A real array of LKJ hyperparameters
   * @param delta A real array of concentration paramters
   * @param shape A vector of shape parameters
   * @param t An integer indicating the number of group-specific terms
   * @param p An integer array with the number variables on the LHS of each |
   * @return target()
   */
  real decov_lp(vector z_b, vector z_T, vector rho, vector zeta, vector tau,
                real[] regularization, real[] delta, vector shape,
                int t, int[] p) {
    int pos_reg = 1;
    int pos_rho = 1;
    target += normal_lpdf(z_b | 0, 1);
    target += normal_lpdf(z_T | 0, 1);
    for (i in 1:t) if (p[i] > 1) {
      vector[p[i] - 1] shape1;
      vector[p[i] - 1] shape2;
      real nu = regularization[pos_reg] + 0.5 * (p[i] - 2);
      pos_reg += 1;
      shape1[1] = nu;
      shape2[1] = nu;
      for (j in 2:(p[i]-1)) {
        nu -= 0.5;
        shape1[j] = 0.5 * j;
        shape2[j] = nu;
      }
      target += beta_lpdf(rho[pos_rho:(pos_rho + p[i] - 2)] | shape1, shape2);
      pos_rho += p[i] - 1;
    }
    target += gamma_lpdf(zeta | delta, 1);
    target += gamma_lpdf(tau  | shape, 1);
    return target();
  }
  
  /**
   * Hierarchical shrinkage parameterization
   *
   * @param z_beta A vector of primitive coefficients
   * @param global A real array of positive numbers
   * @param local A vector array of positive numbers
   * @param global_prior_scale A positive real number
   * @param error_scale 1 or sigma in the Gaussian case
   * @param c2 A positive real number
   * @return A vector of coefficientes
   */
  vector hs_prior(vector z_beta, real[] global, vector[] local, 
                  real global_prior_scale, real error_scale, real c2) {
    int K = rows(z_beta);
    vector[K] lambda = local[1] .* sqrt(local[2]);
    real tau = global[1] * sqrt(global[2]) * global_prior_scale * error_scale;
    vector[K] lambda2 = square(lambda);
    vector[K] lambda_tilde = sqrt( c2 * lambda2 ./ (c2 + square(tau) * lambda2) );
    return z_beta .* lambda_tilde * tau;
  }

  /** 
   * Hierarchical shrinkage plus parameterization
   *
   * @param z_beta A vector of primitive coefficients
   * @param global A real array of positive numbers
   * @param local A vector array of positive numbers
   * @param global_prior_scale A positive real number
   * @param error_scale 1 or sigma in the Gaussian case
   * @param c2 A positive real number
   * @return A vector of coefficientes
   */
  vector hsplus_prior(vector z_beta, real[] global, vector[] local, 
                      real global_prior_scale, real error_scale, real c2) {
    int K = rows(z_beta);
    vector[K] lambda = local[1] .* sqrt(local[2]);
    vector[K] eta = local[3] .* sqrt(local[4]);
    real tau = global[1] * sqrt(global[2]) * global_prior_scale * error_scale;
    vector[K] lambda_eta2 = square(lambda .* eta);
    vector[K] lambda_tilde = sqrt( c2 * lambda_eta2 ./ 
                                 ( c2 + square(tau) * lambda_eta2) );
    return z_beta .* lambda_tilde * tau;
  }
  
  /** 
   * Cornish-Fisher expansion for standard normal to Student t
   *
   * See result 26.7.5 of
   * http://people.math.sfu.ca/~cbm/aands/page_949.htm
   *
   * @param z A scalar distributed standard normal
   * @param df A scalar degrees of freedom
   * @return An (approximate) Student t variate with df degrees of freedom
   */
  real CFt(real z, real df) {
    real z2 = square(z);
    real z3 = z2 * z;
    real z5 = z2 * z3;
    real z7 = z2 * z5;
    real z9 = z2 * z7;
    real df2 = square(df);
    real df3 = df2 * df;
    real df4 = df2 * df2;
    return z + (z3 + z) / (4 * df) + (5 * z5 + 16 * z3 + 3 * z) / (96 * df2)
           + (3 * z7 + 19 * z5 + 17 * z3 - 15 * z) / (384 * df3)
           + (79 * z9 + 776 * z7 + 1482 * z5 - 1920 * z3 - 945 * z) / (92160 * df4);
  }

  /** 
   * Return two-dimensional array of group membership
   *
   * @param N An integer indicating the number of observations
   * @param t An integer indicating the number of grouping variables
   * @param v An integer array with the indices of group membership
   * @return An two-dimensional integer array of group membership
   */
  int[,] make_V(int N, int t, int[] v) {
    int V[t,N];
    int pos = 1;
    if (t > 0) for (j in 1:N) for (i in 1:t) {
      V[i,j] = v[pos] + 1;
      pos += 1;
    }
    return V;
  }

  /**
   * Calculate lower bound on intercept
   *
   * @param family Integer family code
   *   1 = gaussian
   *   2 = gamma
   *   3 = inv-gaussian
   *   4 = beta
   *   5 = binomial
   *   6 = poisson
   *   7 = neg-binom
   *   8 = poisson w/ gamma noise (not currently used but in count.stan)
   * @param link Integer link code
   * @return real lower bound
   */
  real make_lower(int family, int link) {
    if (family == 1) return negative_infinity(); // Gaussian
    if (family <= 3) { // Gamma or inverse Gaussian
      if (link == 2) return negative_infinity(); // log
      return 0;
    }
    return negative_infinity();
  }

  /**
   * Calculate upper bound on intercept
   *
   * @param family Integer family code (see make_lower above for codes)
   * @param link Integer link code
   * @return real upper bound
   */
  real make_upper(int family, int link) {
    if (family == 4 && link == 5) return 0;
    return positive_infinity();
  }
  /** 
   * Apply inverse link function to linear predictor
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_gauss(vector eta, int link) {
    if (link == 1)      return eta;
    else if (link == 2) return exp(eta); 
    else if (link == 3) return inv(eta);
    else reject("Invalid link");
    return eta; // never reached
  }

  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y A vector corresponding to the outcome variable.
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_gauss(vector y, vector eta, real sigma, int link) {
    return -0.5 * log(6.283185307179586232 * sigma) - 
            0.5 * square((y - linkinv_gauss(eta, link)) / sigma);
  }

  /** 
  * Apply inverse link function to linear predictor
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_gamma(vector eta, int link) {
    if (link == 1)      return eta;
    else if (link == 2) return exp(eta);
    else if (link == 3) return inv(eta);
    else reject("Invalid link");
    return eta; // never reached
  }

  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y A vector corresponding to the outcome variable.
  * @param eta A vector of linear predictors
  * @param shape A real number for the shape parameter
  * @param link An integer indicating the link function
  * @param sum_log_y A scalar equal to the sum of log(y)
  * @return A scalar log-likelihood
  */
  real GammaReg(vector y, vector eta, real shape, 
                int link, real sum_log_y) {
    real ret = rows(y) * (shape * log(shape) - lgamma(shape)) +
               (shape - 1) * sum_log_y;
    if (link == 2)      // link is log
      ret -= shape * sum(eta) + shape * sum(y ./ exp(eta));
    else if (link == 1) // link is identity
      ret -= shape * sum(log(eta)) + shape * sum(y ./ eta);
    else if (link == 3) // link is inverse
      ret += shape * sum(log(eta)) - shape * dot_product(eta, y);
    else reject("Invalid link");
    return ret;
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y A vector corresponding to the outcome variable.
  * @param shape A real number for the shape parameter
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_gamma(vector y, vector eta, real shape, int link) {
    int N = rows(eta);
    vector[N] ll;
    if (link == 3) { // link = inverse
      for (n in 1:N) {
        ll[n] = gamma_lpdf(y[n] | shape, shape * eta[n]);
      }
    }
    else if (link == 2) { // link = log
      for (n in 1:N) {
        ll[n] = gamma_lpdf(y[n] | shape, shape / exp(eta[n]));
      }
    }
    else if (link == 1) { // link = identity
      for (n in 1:N) {
        ll[n] = gamma_lpdf(y[n] | shape, shape / eta[n]);
      }
    }
    else reject("Invalid link");
    return ll;
  }

  /** 
  * Apply inverse link function to linear predictor
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_inv_gaussian(vector eta, int link) {
    if (link == 1)      return eta;
    else if (link == 2) return exp(eta);
    else if (link == 3) return inv(eta);
    else if (link == 4) return inv_sqrt(eta);
    else reject("Invalid link");
    return eta; // never reached
  }

  /** 
  * inverse Gaussian log-PDF
  *
  * @param y The vector of outcomes
  * @param mu The vector of conditional means
  * @param lambda A positive scalar dispersion parameter
  * @param sum_log_y A scalar equal to the sum of log(y)
  * @param sqrt_y A vector equal to sqrt(y)
  * @return A scalar
  */
  real inv_gaussian(vector y, vector mu, real lambda, 
                    real sum_log_y, vector sqrt_y) {
    return 0.5 * rows(y) * log(lambda / 6.283185307179586232) - 
      1.5 * sum_log_y - 
      0.5 * lambda * dot_self( (y - mu) ./ (mu .* sqrt_y) );
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y A vector corresponding to the outcome variable.
  * @param eta The linear predictors
  * @param lamba A positive scalar dispersion parameter
  * @param link An integer indicating the link function
  * @param log_y A precalculated vector of the log of y
  * @param sqrt_y A precalculated vector of the square root of y
  * @return A vector of log-likelihoods
  */
  vector pw_inv_gaussian(vector y, vector eta, real lambda, 
                         int link, vector log_y, vector sqrt_y) {
    vector[rows(y)] mu = linkinv_inv_gaussian(eta, link); // link checked
    return -0.5 * lambda * square( (y - mu) ./ (mu .* sqrt_y) ) +
            0.5 * log(lambda / 6.283185307179586232) - 1.5 * log_y;
  }
  
  /** 
  * PRNG for the inverse Gaussian distribution
  *
  * Algorithm from wikipedia 
  *
  * @param mu The expectation
  * @param lambda The dispersion
  * @return A draw from the inverse Gaussian distribution
  */
  real inv_gaussian_rng(real mu, real lambda) {
    real mu2 = square(mu);
    real z = uniform_rng(0,1);
    real y = square(normal_rng(0,1));
    real x = mu + ( mu2 * y - mu * sqrt(4 * mu * lambda * y + mu2 * square(y)) )
           / (2 * lambda);
    if (z <= (mu / (mu + x))) return x;
    else return mu2 / x;
  }
  
  /** 
  * Apply inverse link function to linear predictor for beta models
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_beta(vector eta, int link) {
    if (link == 1) return inv_logit(eta);  // logit
    else if (link == 2) return Phi(eta);   // probit
    else if (link == 3) return inv_cloglog(eta);  // cloglog
    else if (link == 4) return 0.5 + atan(eta) / pi(); // cauchy
    else if (link == 5) return exp(eta); // log 
    else if (link == 6) return 1 - inv_cloglog(-eta); // loglog
    else reject("invalid link");
    return eta; // never reached
  }
  
  /** 
  * Apply inverse link function to linear predictor for dispersion for beta models
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_beta_z(vector eta, int link) {
    if (link == 1) return exp(eta);         // log
    else if (link == 2) return eta;         // identity
    else if (link == 3) return square(eta); // sqrt
    else reject("Invalid link");
    return eta; // never reached
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector for beta models
  *
  * @param y The vector of outcomes
  * @param eta The linear predictors
  * @param dispersion Positive dispersion parameter
  * @param link An integer indicating the link function
  * @return A vector of log-likelihoods
  */
  vector pw_beta(vector y, vector eta, real dispersion, int link) {
    vector[rows(y)] ll;
    vector[rows(y)] mu = linkinv_beta(eta, link); // link checked
    for (n in 1:rows(y)) {
      ll[n] = beta_lpdf(y[n] | mu[n] * dispersion, (1 - mu[n]) * dispersion);
    }
    return ll;
  }

  /** 
  * Pointwise (pw) log-likelihood vector for beta models with z variables
  *
  * @param y The vector of outcomes
  * @param eta The linear predictors (for y)
  * @param eta_z The linear predictors (for dispersion)
  * @param link An integer indicating the link function passed to linkinv_beta
  * @param link_phi An integer indicating the link function passed to linkinv_beta_z
  * @return A vector of log-likelihoods
  */
  vector pw_beta_z(vector y, vector eta, vector eta_z, int link, int link_phi) {
    vector[rows(y)] ll;
    vector[rows(y)] mu = linkinv_beta(eta, link); // link checked
    vector[rows(y)] mu_z = linkinv_beta_z(eta_z, link_phi); // link checked
    for (n in 1:rows(y)) {
      ll[n] = beta_lpdf(y[n] | mu[n] * mu_z[n], (1-mu[n]) * mu_z[n]);
    }
    return ll;
  }
	vector SS_asymp(vector input, matrix Phi_);
	vector SS_asympOff(vector input, matrix Phi_);
	vector SS_asympOrig(vector input, matrix Phi_);
	vector SS_biexp(vector input, matrix Phi_);
	vector SS_fol(vector Dose, vector input, matrix Phi_);
	vector SS_fpl(vector input, matrix Phi_);
	vector SS_gompertz(vector x, matrix Phi_);
	vector SS_logis(vector input, matrix Phi_);
	vector SS_micmen(vector input, matrix Phi_);
	vector SS_weibull(vector x, matrix Phi_);

	/* These functions (without the underscores) are all documented in R
		 See also Appendix C of Pinheiro and Bates
	https://books.google.com/books?id=3TVDAAAAQBAJ&lpg=PR3&dq=Pinheiro%20and%20Bates&pg=PA511#v=onepage&q&f=false
		These functions may be numerically unstable
	*/

	vector SS_asymp(vector input, matrix Phi_) {
		// Phi_[,1] = Asym, Phi_[,2] = R0, Phi_[,3] = lrc
		if (rows(Phi_) > 1) {
			vector[rows(Phi_)] Asym = Phi_[,1];
			return Asym  + (Phi_[,2] - Asym) .* exp(-exp(Phi_[,3]) .* input);
		}
		else {
			real Asym = Phi_[1,1];
			return Asym  + (Phi_[1,2] - Asym) * exp(-exp(Phi_[1,3]) * input);
		}
	}

	vector SS_asympOff(vector input, matrix Phi_) {
		// Phi_[,1] = Asym, Phi_[,2] = lrc, Phi_[,3] = c0
		if (rows(Phi_) > 1) 
			return Phi_[ ,1] .* (1 - exp(-exp(Phi_[ ,2]) .* (input - Phi_[ ,3])));
		else 
			return Phi_[1,1]  * (1 - exp(-exp(Phi_[1,2])  * (input - Phi_[1,3])));
	}

	vector SS_asympOrig(vector input, matrix Phi_) {
		// Phi_[,1] = Asym, Phi_[,2] = lrc
		if (rows(Phi_) > 1)
			return Phi_[ ,1] .* (1 - exp(-exp(Phi_[ ,2]) .* input));
		else 
			return Phi_[1,1]  * (1 - exp(-exp(Phi_[1,2])  * input));
	}

	vector SS_biexp(vector input, matrix Phi_) {
		// Phi_[,1] = A1, Phi_[,2] = lrc1, Phi_[,3] = A2, Phi_[,4] = lrc2
		if (rows(Phi_) > 1)
			return Phi_[ ,1] .* exp(-exp(Phi_[ ,2]) .* input) + 
						 Phi_[ ,3] .* exp(-exp(Phi_[ ,4]) .* input);
		else          
			return Phi_[1,1]  * exp(-exp(Phi_[1,2])  * input) + 
						 Phi_[1,3]  * exp(-exp(Phi_[1,4])  * input);
	}

	vector SS_fol(vector Dose, vector input, matrix Phi_) {
		// Phi_[,1] = lKe, Phi_[,2] = lKa, Phi_[,3] = lCl
		int Phi__rows = rows(Phi_);
		if (Phi__rows > 1) {
			vector[Phi__rows] lKe = Phi_[,1];
			vector[Phi__rows] lKa = Phi_[,2];
			vector[Phi__rows] exp_lKe = exp(lKe);
			vector[Phi__rows] exp_lKa = exp(lKa);
			return Dose .* exp(lKe + lKa - Phi_[,3]) .* 
				(exp(-exp_lKe .* input) - exp(-exp_lKa .* input)) ./ (exp_lKa - exp_lKe);
		}
		else {
			real lKe = Phi_[1,1];
			real lKa = Phi_[1,2];
			real exp_lKe = exp(lKe);
			real exp_lKa = exp(lKa);
			return Dose * exp(lKe + lKa - Phi_[1,3]) .* 
				(exp(-exp_lKe * input) - exp(-exp_lKa * input)) / (exp_lKa - exp_lKe);
		}
	}
		
	vector SS_fpl(vector input, matrix Phi_) {
		// Phi_[,1] = A, Phi_[,2] = B, Phi_[,3] = xmid, Phi_[,4] = scal
		// input is generally data so cannot switch signs
		if (rows(Phi_) > 1) {
			vector[rows(Phi_)] A = Phi_[,1];
			return A + (Phi_[,2] - A) ./ (1 + exp((Phi_[,3] - input) ./ exp(Phi_[,4])));
		}
		else {
			real A = Phi_[1,1];
			return A + rep_vector(Phi_[1,2] - A, rows(input)) 
				./ (1 + exp((Phi_[1,3] - input) / exp(Phi_[1,4])));
		}
	}

	vector SS_gompertz(vector x, matrix Phi_) {
		// Phi_[,1] = Asym, Phi_[,2] = b2, Phi_[,3] = b3
		vector[rows(x)] out;
		if (rows(Phi_) > 1) for (i in 1:rows(x))
			out[i] = Phi_[i,1] * exp(-Phi_[i,2] * Phi_[i,3] ^ x[i]);
		else {
			real Asym = Phi_[1,1];
			real b2 = Phi_[1,2];
			real b3 = Phi_[1,3];
			for (i in 1:rows(x)) out[i] = Asym * exp(-b2 * b3 ^ x[i]);
		}
		return out;
	}

	vector SS_logis(vector input, matrix Phi_) {
		// Phi_[,1] = Asym, Phi_[,2] = xmid, Phi_[,3] = scal
		// input is typically data so cannot switch signs of everything
		if (rows(Phi_) > 1)
			return Phi_[,1] ./ (1 + exp( (Phi_[,2] - input) ./ exp(Phi_[,3])));
		else 
			return rep_vector(Phi_[1,1], rows(input)) ./ 
				(1 + exp( (Phi_[1,2] - input) / exp(Phi_[1,3])));
	}

	vector SS_micmen(vector input, matrix Phi_) {
		// Phi_[,1] = Vm, Phi_[,2] = K
		if (rows(Phi_) > 1)
			return Phi_[ ,1] .* input ./ (Phi_[ ,2] + input);
		else
			return Phi_[1,1]  * input ./ (Phi_[1,2] + input);
	}

	vector SS_weibull(vector x, matrix Phi_) {
		// Phi_[,1] = Asym, Phi_[,2] = Drop, Phi_[,3] = lrc, Phi_[,4] = pwr
		vector[rows(x)] out;
		if (rows(Phi_) > 1) for (i in 1:rows(x))
			out[i] = Phi_[i,1] - Phi_[i,2] * exp(-exp(Phi_[i,3]) * x[i] ^ Phi_[i,4]);
		else {
			real Asym = Phi_[1,1];
			real Drop = Phi_[1,2];
			real lrc = Phi_[1,3];
			real pwr = Phi_[1,4];
			for (i in 1:rows(x)) 
				out[i] = Asym - Drop * exp(-exp(lrc) * x[i] ^ pwr);
		}
		return out;
	}

	matrix reshape_vec(vector x, int Rows, int Cols) {
		matrix[Rows, Cols] out;
		int pos = 1;
		if (rows(x) != Rows * Cols) reject("x is the wrong length");
		for (c in 1:Cols) for (r in 1:Rows) {
			out[r,c] = x[pos];
			pos += 1;
		}
		return out;
	}