#' @noRd
#' @keywords internal
#' @importFrom stats Gamma cor dlogis dnorm gaussian is.empty.model
#' lm median model.matrix model.offset model.response model.weights
#' qbeta qexp rgamma rnorm rt runif sd uniroot
#' @importFrom methods as is
#' @importFrom utils shortPathName

# This function is a slightly modified version of rstanarm::stan_glm.fit
# All credit to https://github.com/stan-dev/stan
#
# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
#
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

cmdstan_glm.fit <-
  function(x, y,
           weights = rep(1, NROW(y)),
           offset = rep(0, NROW(y)),
           family = gaussian(),
           algorithm = c("sampling", "meanfield", "fullrank"),
           out_dir = NULL,
           ...,
           prior = default_prior_coef(family),
           prior_intercept = default_prior_intercept(family),
           prior_aux = exponential(autoscale = TRUE),
           prior_smooth = exponential(autoscale = FALSE),
           prior_PD = FALSE,
           mean_PPD = !prior_PD,
           sparse = FALSE) {

  algorithm <- match.arg(algorithm)
  family <- validate_family(family)
  supported_families <- c("gaussian", "Gamma")
  fam <- which(pmatch(supported_families, family$family, nomatch = 0L) == 1L)
  if (!length(fam)) {
    supported_families_err <- supported_families
    stop("'family' must be one of ", paste(supported_families_err, collapse = ", "))
  }

  supported_links <- supported_glm_links(supported_families[fam])
  link <- which(supported_links == family$link)
  if (!length(link))
    stop("'link' must be one of ", paste(supported_links, collapse = ", "))

  y <- validate_glm_outcome_support(y, family)
  trials <- NULL

  # useless assignments to pass R CMD check
  has_intercept <-
    prior_df <- prior_df_for_intercept <- prior_df_for_aux <- prior_df_for_smooth <-
    prior_dist <- prior_dist_for_intercept <- prior_dist_for_aux <- prior_dist_for_smooth <-
    prior_mean <- prior_mean_for_intercept <- prior_mean_for_aux <- prior_mean_for_smooth <-
    prior_scale <- prior_scale_for_intercept <- prior_scale_for_aux <- prior_scale_for_smooth <-
    prior_autoscale <- prior_autoscale_for_intercept <- prior_autoscale_for_aux <-
    prior_autoscale_for_smooth <- global_prior_scale <- global_prior_df <- slab_df <-
    slab_scale <- xtemp <- xbar <- NULL

  if (is.list(x)) {
    x_stuff <- center_x(x[[1]], sparse)
    smooth_map <- unlist(lapply(1:(length(x) - 1L), FUN = function(j) {
      rep(j, NCOL(x[[j + 1L]]))
    }))
    S <- do.call(cbind, x[-1L])
  }
  else {
    x_stuff <- center_x(x, sparse)
    S <- matrix(NA_real_, nrow = nrow(x), ncol = 0L)
    smooth_map <- integer()
  }
  for (i in names(x_stuff)) # xtemp, xbar, has_intercept
    assign(i, x_stuff[[i]])
  nvars <- ncol(xtemp)

  ok_dists <- nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus",
                    "laplace", "lasso", "product_normal")
  ok_intercept_dists <- ok_dists[1:3]
  ok_aux_dists <- c(ok_dists[1:3], exponential = "exponential")

  # prior distributions
  prior_stuff <- handle_glm_prior(
    prior,
    nvars,
    link = family$link,
    default_scale = 2.5,
    ok_dists = ok_dists
  )
  # prior_{dist, mean, scale, df, dist_name, autoscale},
  # global_prior_df, global_prior_scale, slab_df, slab_scale
  for (i in names(prior_stuff))
    assign(i, prior_stuff[[i]])

  if (isTRUE(is.list(prior_intercept)) &&
      isTRUE(prior_intercept$default)) {
    m_y <- 0
    if (family$family == "gaussian" && family$link == "identity") {
      if (!is.null(y)) m_y <- mean(y) # y can be NULL if prior_PD=TRUE
    }
    prior_intercept$location <- m_y
  }
  prior_intercept_stuff <- handle_glm_prior(
    prior_intercept,
    nvars = 1,
    default_scale = 2.5,
    link = family$link,
    ok_dists = ok_intercept_dists
  )
  # prior_{dist, mean, scale, df, dist_name, autoscale}_for_intercept
  names(prior_intercept_stuff) <- paste0(names(prior_intercept_stuff), "_for_intercept")
  for (i in names(prior_intercept_stuff))
    assign(i, prior_intercept_stuff[[i]])

  prior_aux_stuff <-
    handle_glm_prior(
      prior_aux,
      nvars = 1,
      default_scale = 1,
      link = NULL, # don't need to adjust scale based on logit vs probit
      ok_dists = ok_aux_dists
    )
  # prior_{dist, mean, scale, df, dist_name, autoscale}_for_aux
  names(prior_aux_stuff) <- paste0(names(prior_aux_stuff), "_for_aux")
  if (is.null(prior_aux)) {
    if (prior_PD)
      stop("'prior_aux' cannot be NULL if 'prior_PD' is TRUE.")
    prior_aux_stuff$prior_scale_for_aux <- Inf
  }
  for (i in names(prior_aux_stuff))
    assign(i, prior_aux_stuff[[i]])

  if (ncol(S) > 0) {   # prior_{dist, mean, scale, df, dist_name, autoscale}_for_smooth
    prior_smooth_stuff <-
      handle_glm_prior(
        prior_smooth,
        nvars = max(smooth_map),
        default_scale = 1,
        link = NULL,
        ok_dists = ok_aux_dists
      )

    names(prior_smooth_stuff) <- paste0(names(prior_smooth_stuff), "_for_smooth")
    if (is.null(prior_smooth)) {
      if (prior_PD)
        stop("'prior_smooth' cannot be NULL if 'prior_PD' is TRUE")
      prior_smooth_stuff$prior_scale_for_smooth <- Inf
    }
    for (i in names(prior_smooth_stuff))
      assign(i, prior_smooth_stuff[[i]])

    prior_scale_for_smooth <- array(prior_scale_for_smooth)
  } else {
    prior_dist_for_smooth <- 0L
    prior_mean_for_smooth <- array(NA_real_, dim = 0)
    prior_scale_for_smooth <- array(NA_real_, dim = 0)
    prior_df_for_smooth <- array(NA_real_, dim = 0)
  }

  famname <- supported_families[fam]
  is_bernoulli <- is.binomial(famname) && all(y %in% 0:1) && is.null(trials)
  is_gaussian <- is.gaussian(famname)
  is_gamma <- is.gamma(famname)
  is_continuous <- is_gaussian || is_gamma # always TRUE

  # require intercept for certain family and link combinations
  if (!has_intercept) {
    linkname <- supported_links[link]
    needs_intercept <- !is_gaussian && linkname == "identity" ||
      is_gamma && linkname == "inverse"
    if (needs_intercept)
      stop("To use this combination of family and link ",
           "the model must have an intercept.")
  }

  # allow prior_PD even if no y variable
  if (is.null(y)) {
    if (!prior_PD) {
      stop("Outcome variable must be specified if 'prior_PD' is not TRUE.")
    } else {
      y <- fake_y_for_prior_PD(N = NROW(x), family = family)
      if (is_gaussian &&
          (prior_autoscale || prior_autoscale_for_intercept || prior_autoscale_for_aux)) {
        message("'y' not specified, will assume sd(y)=1 when calculating scaled prior(s). ")
      }
    }
  }

  if (is_gaussian) {
    ss <- sd(y)
    if (prior_dist > 0L && prior_autoscale)
      prior_scale <- ss * prior_scale
    if (prior_dist_for_intercept > 0L && prior_autoscale_for_intercept)
      prior_scale_for_intercept <-  ss * prior_scale_for_intercept
    if (prior_dist_for_aux > 0L && prior_autoscale_for_aux)
      prior_scale_for_aux <- ss * prior_scale_for_aux
  }
  if (prior_dist > 0L && prior_autoscale) {
    min_prior_scale <- 1e-12
    prior_scale <- pmax(min_prior_scale, prior_scale /
                          apply(xtemp, 2L, FUN = function(x) {
                            num.categories <- length(unique(x))
                            x.scale <- 1
                            if (num.categories == 1) {
                              x.scale <- 1
                            } else {
                              x.scale <- sd(x)
                            }
                            return(x.scale)
                          }))
  }
  prior_scale <-
    as.array(pmin(.Machine$double.xmax, prior_scale))
  prior_scale_for_intercept <-
    min(.Machine$double.xmax, prior_scale_for_intercept)

  if (length(weights) > 0 && all(weights == 1)) weights <- double()
  if (length(offset)  > 0 && all(offset  == 0)) offset  <- double()

  # create entries in the data block of the .stan file
  standata <- nlist(
    N = nrow(xtemp),
    K = ncol(xtemp),
    xbar = as.array(xbar),
    dense_X = !sparse,
    family = stan_family_number(famname),
    link,
    has_weights = length(weights) > 0,
    has_offset = length(offset) > 0,
    has_intercept,
    prior_PD,
    compute_mean_PPD = mean_PPD,
    prior_dist,
    prior_mean,
    prior_scale,
    prior_df,
    prior_dist_for_intercept,
    prior_scale_for_intercept = c(prior_scale_for_intercept),
    prior_mean_for_intercept = c(prior_mean_for_intercept),
    prior_df_for_intercept = c(prior_df_for_intercept),
    global_prior_df, global_prior_scale, slab_df, slab_scale, # for hs priors
    z_dim = 0,  # betareg data
    link_phi = 0,
    betareg_z = array(0, dim = c(nrow(xtemp), 0)),
    has_intercept_z = 0,
    zbar = array(0, dim = c(0)),
    prior_dist_z = 0, prior_mean_z = integer(), prior_scale_z = integer(),
    prior_df_z = integer(), global_prior_scale_z = 0, global_prior_df_z = 0,
    prior_dist_for_intercept_z = 0, prior_mean_for_intercept_z = 0,
    prior_scale_for_intercept_z = 0, prior_df_for_intercept_z = 0,
    prior_df_for_intercept = c(prior_df_for_intercept),
    prior_dist_for_aux = prior_dist_for_aux,
    prior_dist_for_smooth, prior_mean_for_smooth, prior_scale_for_smooth, prior_df_for_smooth,
    slab_df_z = 0, slab_scale_z = 0,
    num_normals = if(prior_dist == 7) as.integer(prior_df) else integer(0),
    num_normals_z = integer(0),
    clogit = 0L, J = 0L, strata = integer()
    # mean,df,scale for aux added below depending on family
  )

  # make a copy of user specification before modifying 'group' (used for keeping
  # track of priors)
  user_covariance <- NULL
  standata$t <- 0L
  standata$p <- integer(0)
  standata$l <- integer(0)
  standata$q <- 0L
  standata$len_theta_L <- 0L
  standata$num_non_zero <- 0L
  standata$w <- double(0)
  standata$v <- integer(0)
  standata$u <- integer(0)
  standata$special_case <- 0L
  standata$shape <- standata$scale <- standata$concentration <-
    standata$regularization <- rep(0, 0)
  standata$len_concentration <- 0L
  standata$len_regularization <- 0L
  standata$SSfun <- 0L
  standata$input <- double()
  standata$Dose <- double()

  if (!is_bernoulli) {
    standata$X <- array(xtemp, dim = c(1L, dim(xtemp)))
    standata$nnz_X <- 0L
    standata$w_X <- double(0)
    standata$v_X <- integer(0)
    standata$u_X <- integer(0)
    standata$y <- y
    standata$weights <- weights
    standata$offset_ <- offset
    standata$K_smooth <- ncol(S)
    standata$S <- S
    standata$smooth_map <- smooth_map
  }

  # call stan() to draw from posterior distribution
  if (is_continuous) {
    standata$ub_y <- Inf
    standata$lb_y <- if (is_gaussian) -Inf else 0
    standata$prior_scale_for_aux <- prior_scale_for_aux %ORifINF% 0
    standata$prior_df_for_aux <- c(prior_df_for_aux)
    standata$prior_mean_for_aux <- c(prior_mean_for_aux)
    standata$len_y <- length(y)
    stan_model <- cmdstanr::cmdstan_model(
      system.file("extdata/stan_files/from_rstanarm/continuous.stan", package = "pstpipeline"),
      include_paths = shortPathName(system.file(
        "extdata/stan_files/from_rstanarm", package = "pstpipeline") ## this will fail if path has spaces!
      )
    )
  }
  else {
    stop(paste(famname, "is not supported in this version of the function."))
  }

  prior_info <- summarize_glm_prior(
    user_prior = prior_stuff,
    user_prior_intercept = prior_intercept_stuff,
    user_prior_aux = prior_aux_stuff,
    user_prior_covariance = user_covariance,
    has_intercept = has_intercept,
    has_predictors = nvars > 0,
    adjusted_prior_scale = prior_scale,
    adjusted_prior_intercept_scale = prior_scale_for_intercept,
    adjusted_prior_aux_scale = prior_scale_for_aux,
    family = family
  )

  pars <- c(if (has_intercept) "alpha",
            "beta",
            if (ncol(S)) "beta_smooth",
            if (is_continuous) "aux",
            if (ncol(S)) "smooth_sd",
            if (standata$len_theta_L) "theta_L",
            if (mean_PPD) "mean_PPD")

  l <- list(...)
  if (algorithm != "sampling") {
    if (is.null(l$iter)) l$iter <- 10000
    if (is.null(l$output_samples)) l$output_samples <- 1000
  }
  else { # clearly nothing is being changed, given here just to show defaults
    if (is.null(l$chains)) l$chains <- 4 # default (explicitly defined here for file naming)
    if (is.null(l$iter_sampling)) l$iter_sampling <- 1000 # default (explicitly defined here for file naming)
    if (is.null(l$adapt_delta)) l$adapt_delta <- 0.95
  }

  cmdstanr::check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
  if (!is.null(out_dir)) {
    out_dir <- file.path(getwd(), out_dir)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  }

  if (algorithm != "sampling") {
    fit <- stan_model$variational(
      data = standata,
      algorithm = algorithm,
      seed = l$seed,
      iter = l$iter,
      refresh = l$refresh,
      output_samples = l$output_samples,
      output_dir = out_dir
    )
  }
  else {
    # MCMC
    fit <- stan_model$sample(
      data = standata,
      seed = l$seed,
      refresh = l$refresh, # default = 100
      chains = l$chains, # default = 4
      iter_warmup = l$iter_warmup, # default = 1000
      iter_sampling = l$iter_sampling, # default = 1000
      adapt_delta = l$adapt_delta, # default = 0.95
      step_size = l$step_size, # default = 1
      max_treedepth = l$max_treedepth, # default = 10
      output_dir = out_dir
    )
  }

  return(fit)

}


# rstanarm::priors internal fns

normal <- function(location = 0, scale = NULL, autoscale = FALSE) {
  validate_parameter_value(scale)
  nlist(dist = "normal", df = NA, location, scale, autoscale)
}
student_t <- function(df = 1, location = 0, scale = NULL, autoscale = FALSE) {
  validate_parameter_value(scale)
  validate_parameter_value(df)
  nlist(dist = "t", df, location, scale, autoscale)
}
cauchy <- function(location = 0, scale = NULL, autoscale = FALSE) {
  student_t(df = 1, location = location, scale = scale, autoscale)
}
hs <- function(df = 1, global_df = 1, global_scale = 0.01,
               slab_df = 4, slab_scale = 2.5) {
  validate_parameter_value(df)
  validate_parameter_value(global_df)
  validate_parameter_value(global_scale)
  validate_parameter_value(slab_df)
  validate_parameter_value(slab_scale)
  nlist(dist = "hs", df, location = 0, scale = 1,
        global_df, global_scale, slab_df, slab_scale)
}
hs_plus <- function(df1 = 1, df2 = 1, global_df = 1, global_scale = 0.01,
                    slab_df = 4, slab_scale = 2.5) {
  validate_parameter_value(df1)
  validate_parameter_value(df2)
  validate_parameter_value(global_df)
  validate_parameter_value(global_scale)
  validate_parameter_value(slab_df)
  validate_parameter_value(slab_scale)
  # scale gets used as a second df hyperparameter
  nlist(dist = "hs_plus", df = df1, location = 0, scale = df2, global_df,
        global_scale, slab_df, slab_scale)
}
laplace <- function(location = 0, scale = NULL, autoscale = FALSE) {
  nlist(dist = "laplace", df = NA, location, scale, autoscale)
}
lasso <- function(df = 1, location = 0, scale = NULL, autoscale = FALSE) {
  nlist(dist = "lasso", df, location, scale, autoscale)
}
product_normal <- function(df = 2, location = 0, scale = 1) {
  validate_parameter_value(df)
  stopifnot(all(df >= 1), all(df == as.integer(df)))
  validate_parameter_value(scale)
  nlist(dist = "product_normal", df, location, scale)
}
exponential <- function(rate = 1, autoscale = FALSE) {
  stopifnot(length(rate) == 1)
  validate_parameter_value(rate)
  nlist(dist = "exponential",
        df = NA, location = NA, scale = 1/rate,
        autoscale)
}
decov <- function(regularization = 1, concentration = 1,
                  shape = 1, scale = 1) {
  validate_parameter_value(regularization)
  validate_parameter_value(concentration)
  validate_parameter_value(shape)
  validate_parameter_value(scale)
  nlist(dist = "decov", regularization, concentration, shape, scale)
}
lkj <- function(regularization = 1, scale = 10, df = 1, autoscale = TRUE) {
  validate_parameter_value(regularization)
  validate_parameter_value(scale)
  validate_parameter_value(df)
  nlist(dist = "lkj", regularization, scale, df, autoscale)
}
dirichlet <- function(concentration = 1) {
  validate_parameter_value(concentration)
  nlist(dist = "dirichlet", concentration)
}
R2 <- function(location = NULL, what = c("mode", "mean", "median", "log")) {
  what <- match.arg(what)
  validate_R2_location(location, what)
  list(dist = "R2", location = location, what = what, df = 0, scale = 0)
}
default_prior_intercept = function(family) {
  # family arg not used, but we can use in the future to do different things
  # based on family if necessary
  out <- normal(0, 2.5, autoscale = TRUE)
  out$location <- NULL # not determined yet
  out$default <- TRUE
  out$version <- utils::packageVersion("rstanarm")
  out
}
default_prior_coef = function(family) {
  # family arg not used, but we can use in the future to do different things
  # based on family if necessary
  out <- normal(0, 2.5, autoscale = TRUE)
  out$default <- TRUE
  out$version <- utils::packageVersion("rstanarm")
  out
}
validate_parameter_value <- function(x) {
  nm <- deparse(substitute(x))
  if (!is.null(x)) {
    if (!is.numeric(x))
      stop(nm, " should be NULL or numeric", call. = FALSE)
    if (any(x <= 0))
      stop(nm, " should be positive", call. = FALSE)
  }
  invisible(TRUE)
}

validate_R2_location <- function(location = NULL, what) {
  stopifnot(is.numeric(location))
  if (length(location) > 1)
    stop(
      "The 'R2' function only accepts a single value for 'location', ",
      "which applies to the prior R^2. ",
      "If you are trying to put different priors on different coefficients ",
      "rather than specify a joint prior via 'R2', you can use stan_glm ",
      "which accepts a wider variety of priors, many of which allow ",
      "specifying arguments as vectors.",
      call. = FALSE
    )

  if (what == "log") {
    if (location >= 0)
      stop("If 'what' is 'log' then location must be negative.", call. = FALSE)
  } else if (what == "mode") {
    if (location <= 0 || location > 1)
      stop("If 'what' is 'mode', location must be in (0,1].",
           call. = FALSE)
  } else { # "mean", "median"
    if (location <= 0 || location >= 1)
      stop("If 'what' is 'mean' or 'median', location must be in (0,1).",
           call. = FALSE)
  }
  invisible(TRUE)
}

make_eta <- function(location, what = c("mode", "mean", "median", "log"), K) {
  stopifnot(length(location) == 1, is.numeric(location))
  stopifnot(is.numeric(K), K == as.integer(K))
  if (K == 0)
    stop("R2 prior is not applicable when there are no covariates.",
         call. = FALSE)
  what <- match.arg(what)
  half_K <- K / 2
  if (what == "mode") {
    stopifnot(location > 0, location <= 1)
    if (K <= 2)
      stop(paste("R2 prior error.",
                 "The mode of the beta distribution does not exist",
                 "with fewer than three predictors.",
                 "Specify 'what' as 'mean', 'median', or 'log' instead."),
           call. = FALSE)
    eta <- (half_K - 1  - location * half_K + location * 2) / location
  } else if (what == "mean") {
    stopifnot(location > 0, location < 1)
    eta <- (half_K - location * half_K) / location
  } else if (what == "median") {
    stopifnot(location > 0, location < 1)
    FUN <- function(eta) qbeta(0.5, half_K, qexp(eta)) - location
    eta <- qexp(uniroot(FUN, interval = 0:1)$root)
  } else { # what == "log"
    stopifnot(location < 0)
    FUN <- function(eta) digamma(half_K) - digamma(half_K + qexp(eta)) - location
    eta <- qexp(uniroot(FUN, interval = 0:1,
                        f.lower = -location,
                        f.upper = -.Machine$double.xmax)$root)
  }

  return(eta)
}

# rstanarm::stan_glm.fit internal fns ----------------------------------------------------------------

validate_family <- function(f) {
  if (is.character(f))
    f <- get(f, mode = "function", envir = parent.frame(2))
  if (is.function(f))
    f <- f()
  if (!is(f, "family"))
    stop("'family' must be a family.", call. = FALSE)

  return(f)
}
center_x <- function(x, sparse) {
  x <- as.matrix(x)
  has_intercept <- if (ncol(x) == 0)
    FALSE else grepl("(Intercept", colnames(x)[1L], fixed = TRUE)

  xtemp <- if (has_intercept) x[, -1L, drop=FALSE] else x
  if (has_intercept && !sparse) {
    xbar <- colMeans(xtemp)
    xtemp <- sweep(xtemp, 2, xbar, FUN = "-")
  }
  else xbar <- rep(0, ncol(xtemp))

  sel <- apply(xtemp, 2L, function(x) !all(x == 1) && length(unique(x)) < 2)
  if (any(sel)) {
    # drop any column of x with < 2 unique values (empty interaction levels)
    # exception is column of 1s isn't dropped
    warning("Dropped empty interaction levels: ",
            paste(colnames(xtemp)[sel], collapse = ", "))
    xtemp <- xtemp[, !sel, drop = FALSE]
    xbar <- xbar[!sel]
  }

  return(nlist(xtemp, xbar, has_intercept))
}
nlist <- function(...) {
  m <- match.call()
  out <- list(...)
  no_names <- is.null(names(out))
  has_name <- if (no_names) FALSE else nzchar(names(out))
  if (all(has_name))
    return(out)
  nms <- as.character(m)[-1L]
  if (no_names) {
    names(out) <- nms
  } else {
    names(out)[!has_name] <- nms[!has_name]
  }

  return(out)
}
is.gaussian <- function(x) x == "gaussian"
is.gamma <- function(x) x == "Gamma"
is.binomial <- function(x) x == "binomial"
`%ORifNULL%` <- function(a, b) {
  if (is.null(a)) b else a
}
`%ORifINF%` <- function(a, b) {
  if (a == Inf) b else a
}
set_prior_scale <- function(scale, default, link) {
  stopifnot(is.numeric(default), is.character(link) || is.null(link))
  if (is.null(scale))
    scale <- default
  if (isTRUE(link == "probit"))
    scale <- scale * dnorm(0) / dlogis(0)

  return(scale)
}
drop_redundant_dims <- function(data) {
  drop_dim <- sapply(data, function(v) is.matrix(v) && NCOL(v) == 1)
  data[, drop_dim] <- lapply(data[, drop_dim, drop=FALSE], drop)
  return(data)
}
validate_data <- function(data, if_missing = NULL) {
  if (missing(data) || is.null(data)) {
    warning("Omitting the 'data' argument is not recommended.")
    return(if_missing)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }

  # drop other classes (e.g. 'tbl_df', 'tbl', 'data.table')
  data <- as.data.frame(data)
  drop_redundant_dims(data)
}
check_constant_vars <- function(mf) {
  mf1 <- mf
  if (NCOL(mf[, 1]) == 2 || all(mf[, 1] %in% c(0, 1))) {
    mf1 <- mf[, -1, drop=FALSE]
  }

  lu1 <- function(x) !all(x == 1) && length(unique(x)) == 1
  nocheck <- c("(weights)", "(offset)", "(Intercept)")
  sel <- !colnames(mf1) %in% nocheck
  is_constant <- apply(mf1[, sel, drop=FALSE], 2, lu1)
  if (any(is_constant)) {
    stop("Constant variable(s) found: ",
         paste(names(is_constant)[is_constant], collapse = ", "),
         call. = FALSE)
  }
  return(mf)
}
handle_glm_prior <- function(prior, nvars, default_scale, link,
                             ok_dists = nlist("normal", student_t = "t",
                                              "cauchy", "hs", "hs_plus",
                                              "laplace", "lasso", "product_normal")) {
  if (!length(prior))
    return(list(prior_dist = 0L, prior_mean = as.array(rep(0, nvars)),
                prior_scale = as.array(rep(1, nvars)),
                prior_df = as.array(rep(1, nvars)), prior_dist_name = NA,
                global_prior_scale = 0, global_prior_df = 0,
                slab_df = 0, slab_scale = 0,
                prior_autoscale = FALSE))

  if (!is.list(prior))
    stop(sQuote(deparse(substitute(prior))), " should be a named list")

  prior_dist_name <- prior$dist
  prior_scale <- prior$scale
  prior_mean <- prior$location
  prior_df <- prior$df
  prior_mean[is.na(prior_mean)] <- 0
  prior_df[is.na(prior_df)] <- 1
  global_prior_scale <- 0
  global_prior_df <- 0
  slab_df <- 0
  slab_scale <- 0
  if (!prior_dist_name %in% unlist(ok_dists)) {
    stop("The prior distribution should be one of ",
         paste(names(ok_dists), collapse = ", "))
  } else if (prior_dist_name %in%
             c("normal", "t", "cauchy", "laplace", "lasso", "product_normal")) {
    if (prior_dist_name == "normal") prior_dist <- 1L
    else if (prior_dist_name == "t") prior_dist <- 2L
    else if (prior_dist_name == "laplace") prior_dist <- 5L
    else if (prior_dist_name == "lasso") prior_dist <- 6L
    else if (prior_dist_name == "product_normal") prior_dist <- 7L
    prior_scale <- set_prior_scale(prior_scale, default = default_scale,
                                   link = link)
  } else if (prior_dist_name %in% c("hs", "hs_plus")) {
    prior_dist <- ifelse(prior_dist_name == "hs", 3L, 4L)
    global_prior_scale <- prior$global_scale
    global_prior_df <- prior$global_df
    slab_df <- prior$slab_df
    slab_scale <- prior$slab_scale
  } else if (prior_dist_name %in% "exponential") {
    prior_dist <- 3L # only used for scale parameters so 3 not a conflict with 3 for hs
  }

  prior_df <- maybe_broadcast(prior_df, nvars)
  prior_df <- as.array(pmin(.Machine$double.xmax, prior_df))
  prior_mean <- maybe_broadcast(prior_mean, nvars)
  prior_mean <- as.array(prior_mean)
  prior_scale <- maybe_broadcast(prior_scale, nvars)

  nlist(prior_dist,
        prior_mean,
        prior_scale,
        prior_df,
        prior_dist_name,
        global_prior_scale,
        global_prior_df,
        slab_df,
        slab_scale,
        prior_autoscale = isTRUE(prior$autoscale))
}
supported_glm_links <- function(famname) {
  switch(
    famname,
    gaussian = c("identity", "log", "inverse"),
    Gamma = c("identity", "log", "inverse"),
    stop("unsupported family")
  )
}
stan_family_number <- function(famname) {
  switch(
    famname,
    "gaussian" = 1L,
    "Gamma" = 2L,
    stop("Family not valid.")
  )
}
validate_glm_outcome_support <- function(y, family) {
  if (is.character(y)) {
    stop("Outcome variable can't be type 'character'.", call. = FALSE)
  }

  if (is.null(y)) {
    return(y)
  }

  .is_count <- function(x) {
    all(x >= 0) && all(abs(x - round(x)) < .Machine$double.eps^0.5)
  }

  fam <- family$family

  if (!is.binomial(fam)) {
    # make sure y has ok dimensions (matrix only allowed for binomial models)
    if (length(dim(y)) > 1) {
      if (NCOL(y) == 1) {
        y <- y[, 1]
      } else {
        stop("Except for binomial models the outcome variable ",
             "should not have multiple columns.",
             call. = FALSE)
      }
    }
    # check that values match support for non-binomial models
    if (is.gaussian(fam)) {
      return(y)
    }
    else if (is.gamma(fam) && any(y <= 0)) {
      stop("All outcome values must be positive for gamma models.",
           call. = FALSE)
    }
  }

  return(y)
}
fake_y_for_prior_PD <- function(N, family) {
  fam <- family$family
  if (is.gaussian(fam)) {
    # if prior autoscaling is on then the value of sd(y) matters
    # generate a fake y so that sd(y) is 1
    fake_y <- as.vector(scale(rnorm(N)))
  }
  else {
    # valid for gamma, inverse gaussian, beta
    fake_y <- runif(N)
  }
  return(fake_y)
}

# Create "prior.info" attribute needed for prior_summary()
summarize_glm_prior <-
  function(user_prior,
           user_prior_intercept,
           user_prior_aux,
           user_prior_covariance,
           has_intercept,
           has_predictors,
           adjusted_prior_scale,
           adjusted_prior_intercept_scale,
           adjusted_prior_aux_scale,
           family) {
    rescaled_coef <-
      user_prior$prior_autoscale &&
      has_predictors &&
      !is.na(user_prior$prior_dist_name) &&
      !all(user_prior$prior_scale == adjusted_prior_scale)
    rescaled_int <-
      user_prior_intercept$prior_autoscale_for_intercept &&
      has_intercept &&
      !is.na(user_prior_intercept$prior_dist_name_for_intercept) &&
      (user_prior_intercept$prior_scale_for_intercept != adjusted_prior_intercept_scale)
    rescaled_aux <- user_prior_aux$prior_autoscale_for_aux &&
      !is.na(user_prior_aux$prior_dist_name_for_aux) &&
      (user_prior_aux$prior_scale_for_aux != adjusted_prior_aux_scale)

    if (has_predictors && user_prior$prior_dist_name %in% "t") {
      if (all(user_prior$prior_df == 1)) {
        user_prior$prior_dist_name <- "cauchy"
      } else {
        user_prior$prior_dist_name <- "student_t"
      }
    }
    if (has_intercept &&
        user_prior_intercept$prior_dist_name_for_intercept %in% "t") {
      if (all(user_prior_intercept$prior_df_for_intercept == 1)) {
        user_prior_intercept$prior_dist_name_for_intercept <- "cauchy"
      } else {
        user_prior_intercept$prior_dist_name_for_intercept <- "student_t"
      }
    }
    if (user_prior_aux$prior_dist_name_for_aux %in% "t") {
      if (all(user_prior_aux$prior_df_for_aux == 1)) {
        user_prior_aux$prior_dist_name_for_aux <- "cauchy"
      } else {
        user_prior_aux$prior_dist_name_for_aux <- "student_t"
      }
    }
    prior_list <- list(
      prior =
        if (!has_predictors) NULL else with(user_prior, list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_coef)
            adjusted_prior_scale else NULL,
          df = if (prior_dist_name %in% c
                   ("student_t", "hs", "hs_plus", "lasso", "product_normal"))
            prior_df else NULL
        )),
      prior_intercept =
        if (!has_intercept) NULL else with(user_prior_intercept, list(
          dist = prior_dist_name_for_intercept,
          location = prior_mean_for_intercept,
          scale = prior_scale_for_intercept,
          adjusted_scale = if (rescaled_int)
            adjusted_prior_intercept_scale else NULL,
          df = if (prior_dist_name_for_intercept %in% "student_t")
            prior_df_for_intercept else NULL
        ))
    )
    if (length(user_prior_covariance))
      prior_list$prior_covariance <- user_prior_covariance

    aux_name <- .rename_aux(family)
    prior_list$prior_aux <- if (is.na(aux_name))
      NULL else with(user_prior_aux, list(
        dist = prior_dist_name_for_aux,
        location = if (!is.na(prior_dist_name_for_aux) &&
                       prior_dist_name_for_aux != "exponential")
          prior_mean_for_aux else NULL,
        scale = if (!is.na(prior_dist_name_for_aux) &&
                    prior_dist_name_for_aux != "exponential")
          prior_scale_for_aux else NULL,
        adjusted_scale = if (rescaled_aux)
          adjusted_prior_aux_scale else NULL,
        df = if (!is.na(prior_dist_name_for_aux) &&
                 prior_dist_name_for_aux %in% "student_t")
          prior_df_for_aux else NULL,
        rate = if (!is.na(prior_dist_name_for_aux) &&
                   prior_dist_name_for_aux %in% "exponential")
          1 / prior_scale_for_aux else NULL,
        aux_name = aux_name
      ))

    return(prior_list)
  }
array1D_check <- function(y) {
  if (length(dim(y)) == 1L) {
    nms <- rownames(y)
    dim(y) <- NULL
    if (!is.null(nms))
      names(y) <- nms
  }
  return(y)
}
# rename aux parameter based on family
.rename_aux <- function(family) {
  fam <- family$family
  if (is.gaussian(fam)) "sigma" else
    if (is.gamma(fam)) "shape"
}
.sample_indices <- function(wts, n_draws) {
  ## Stratified resampling
  ##   Kitagawa, G., Monte Carlo Filter and Smoother for Non-Gaussian
  ##   Nonlinear State Space Models, Journal of Computational and
  ##   Graphical Statistics, 5(1):1-25, 1996.
  K <- length(wts)
  w <- n_draws * wts # expected number of draws from each model
  idx <- rep(NA, n_draws)

  c <- 0
  j <- 0

  for (k in 1:K) {
    c <- c + w[k]
    if (c >= 1) {
      a <- floor(c)
      c <- c - a
      idx[j + 1:a] <- k
      j <- j + a
    }
    if (j < n_draws && c >= stats::runif(1)) {
      c <- c - 1
      j <- j + 1
      idx[j] <- k
    }
  }
  return(idx)
}
maybe_broadcast <- function(x, n) {
  if (!length(x)) {
    rep(0, times = n)
  } else if (length(x) == 1L) {
    rep(x, times = n)
  } else {
    x
  }
}
