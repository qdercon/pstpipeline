#' Helper function to take a subsample of our parsed list (for demonstration
#' purposes)
#'
#' \code{take_subsample} takes a sample of size n to enable model fitting on a
#' smaller sample.
#'
#' @param parsed_list A list outputted from [import_multiple()].
#' @param n_ppts Sample size to take.
#'
#' @returns A named \code{list}.
#'
#' @examples
#' data(example_data)
#' subsamp <- take_subsample(example_data$nd, 10) # sample of 10 participants
#' @export

take_subsample <- function(parsed_list,
                           n_ppts) {

  if (is.null(parsed_list[["ppt_info"]])) {
    stop(
      strwrap("Could not find a list of participant info to take a sample of
              IDs. Perhaps the list is split?", prefix = " ", initial = "")
    )
  }

  subjID <- NULL # to appease R CMD check
  ids <- sample(unique(parsed_list[["ppt_info"]][["subjID"]]), size = n_ppts)
  subsample <- list()
  elements <- names(parsed_list)

  for (el in elements) {
    subsample[[el]] <- parsed_list[[el]] |>
      dplyr::filter(subjID %in% ids)
  }

  subsample
}

#' Example probabilistic selection task data
#'
#' An example dataset with data from ten individuals from each group.
#'
#' @docType data
#' @keywords example_data
#' @name example_data
#' @usage data(example_data)
#' @format A list with two elements: \code{nd} has ten non-distanced
#' individuals, and \code{dis} has ten distanced individuals. Each sub-list is a
#' list with four elements: (\code{ppt_info}, \code{training}, \code{test}, and
#' \code{gillan_questions}).
NULL

#' Compute a standard error of the mean
#'
#' \code{std} computes the standard error of a mean.
#'
#' @param x A numerical vector of values.
#'
#' @returns A numeric value.
#' @noRd

std <- function(x) sd(x, na.rm = TRUE) / sqrt(length(x))

#' Parse ELBO at a target iteration from cmdstanr variational latent dynamics
#'
#' @param latent_path Path to latent dynamics CSV file.
#' @param target_iter Iteration for which ELBO should be extracted.
#'
#' @returns A named list with parsed ELBO and debug context.
#' @noRd

extract_elbo_from_file <- function(latent_path,
                                   target_iter = 100) {
  if (length(latent_path) == 0 || is.na(latent_path) || latent_path == "") {
    return(
      list(
        target_elbo = NA_real_,
        parsed_rows = 0,
        lines = character(0),
        file = latent_path
      )
    )
  }

  if (!file.exists(latent_path)) {
    return(
      list(
        target_elbo = NA_real_,
        parsed_rows = 0,
        lines = character(0),
        file = latent_path
      )
    )
  }

  lines <- readLines(latent_path, warn = FALSE)
  dat <- tryCatch(
    utils::read.csv(
      latent_path,
      header = FALSE,
      comment.char = "#",
      stringsAsFactors = FALSE
    ),
    error = function(e) data.frame()
  )

  if (!is.data.frame(dat) || nrow(dat) == 0 || ncol(dat) < 2) {
    return(
      list(
        target_elbo = NA_real_,
        parsed_rows = 0,
        lines = lines,
        file = latent_path
      )
    )
  }

  iter_vals <- suppressWarnings(as.integer(dat[[1]]))
  if (ncol(dat) >= 3) {
    elbo_vals <- suppressWarnings(as.numeric(dat[[3]]))
  } else {
    elbo_vals <- suppressWarnings(as.numeric(dat[[2]]))
  }

  ok <- !is.na(iter_vals) & !is.na(elbo_vals)

  if (!any(ok)) {
    return(
      list(
        target_elbo = NA_real_,
        parsed_rows = 0,
        lines = lines,
        file = latent_path
      )
    )
  }

  parsed_df <- data.frame(
    iter = iter_vals[ok],
    elbo = elbo_vals[ok],
    stringsAsFactors = FALSE
  )

  target_idx <- which(parsed_df$iter == target_iter)
  target_elbo <- if (length(target_idx) > 0) {
    parsed_df$elbo[target_idx[length(target_idx)]]
  } else {
    NA_real_
  }

  list(
    target_elbo = target_elbo,
    parsed_rows = nrow(parsed_df),
    lines = lines,
    file = latent_path
  )
}

#' Compute a single highest posterior density interval (HDI)
#'
#' \code{single_hdi} computes the highest density interval from a sample of
#' representative values, estimated as shortest credible interval.
#'
#' @param vals A vector of representative values from a probability distribution
#' (e.g., MCMC samples).
#' @param cred A scalar between 0 and 1, indicating the mass within the
#' credible interval that is to be estimated.
#'
#' @returns A vector with the lower and upper HDI.
#' @noRd

# Adapted from hBayesDM::HDIofMCMC, in turn based on John Kruschke's code.

single_hdi <- function(vals,
                       cred) {
  sampleVec <- as.vector(t(vals))
  sortedPts <- sort(sampleVec)
  ciIdxInc <- floor(cred * length(sortedPts))
  nCIs <- length(sortedPts) - ciIdxInc
  ciWidth <- rep(0, nCIs)
  for (i in 1:nCIs) {
    ciWidth[i] <- sortedPts[i + ciIdxInc] - sortedPts[i]
  }
  HDImin <- sortedPts[which.min(ciWidth)]
  HDImax <- sortedPts[which.min(ciWidth) + ciIdxInc]
  HDIlim <- c(HDImin, HDImax)
  as.vector(t(HDIlim))
}

#' Compute quantiles of a probability distrbution based on highest density
#' intervals (HDIs)
#'
#' \code{quantile_hdi} computes any number of highest density quantiles from a
#' sample of representative values, estimated as shortest credible intervals. If
#' \code{quantile} contains 0, 0.5, or 1, will return the minimum, median, and
#' maximum respectively.
#'
#' @param var A vector of representative values from a probability distribution
#' (e.g., MCMC samples).
#' @param quantile A vector of quantiles to return.
#' @param moment A string indicating the first moment to return ("mean" or
#' "median") if \code{quantile} contains 0.5.
#' @param ... Internal arguments.
#'
#' @returns A sorted vector with all specified quantiles.
#'
#' @examples
#' p_density <- rnorm(100, 2, 0.5)
#' quantile_hdi(p_density, c(0, 0.025, 0.5, 0.0975, 1))
#' @export

quantile_hdi <- function(var,
                         quantile,
                         moment = "mean",
                         ...) {

  l <- list(...)
  if (is.null(l$transform)) l$transform <- FALSE
  cent <- match.arg(moment, c("mean", "median"), several.ok = FALSE)
  if (l$transform) {
    # fixes issues calling fn from within ggplot with exponentiated coefficients
    var <- log(var / 100 + 1)
  }
  returns <- vector(mode = "numeric")

  for (q in quantile) {
    if (q == 0.5) {
      returns <- cbind(returns, get(cent)(var))
    } else if (q == 0) {
      returns <- cbind(returns, min(var))
    } else if (q == 1) {
      returns <- cbind(returns, max(var))
    } else if (q < 0.5) {
      cred_mass <- 1 - 2 * q
      HDI_lower <- single_hdi(vals = var, cred = cred_mass)[1]
      returns <- cbind(returns, HDI_lower)
    } else {
      cred_mass <- 2 * q - 1
      HDI_upper <- single_hdi(vals = var, cred = cred_mass)[2]
      returns <- cbind(returns, HDI_upper)
    }
  }

  ret <- sort(returns)
  if (l$transform) {
    ret <- (exp(ret) - 1) * 100
  }
  names(ret) <- sapply(
    seq_along(quantile), FUN = function(x) paste0(quantile[x] * 100, "%")
  )

  ret
}

#' Define GLM family based on parameter name
#'
#' \code{family_ch} returns \code{Gamma(link = "log")} if the parameter is a
#' learning rate (\code{alpha}) or decay factor (\code{gamma}), or
#' \code{gaussian()} for any other parameter name.
#'
#' @param param Parameter name.
#'
#' @returns A [stats::family()].
#' @noRd

family_ch <- function(param) {
  if (grepl("alpha|gamma", param)) Gamma(link = "log")
  else gaussian()
}

#' Identify individual-level learning-model parameters among Stan variable
#' names
#'
#' \code{is_learning_model_var} matches the (base or bracket-indexed) variable
#' names used throughout the package (e.g., in \code{clean_summary}), excluding
#' group-level hyperparameters and generated quantities such as \code{y_pred}
#' and \code{log_lik}.
#'
#' @param x Character vector of variable names.
#'
#' @returns A logical vector, the same length as \code{x}.
#' @noRd

is_learning_model_var <- function(x) {
  grepl("alpha|beta|rho|pers|kappa|w|gamma", x) &
    !grepl("_pr|_s|mu|sigma|_diff|_i|pointwise", x)
}

#' Clean up summary output from cmdstanr
#'
#' \code{clean_summary} returns a long format [tibble::tibble()] with individual
#' parameters from learning modelsl.
#'
#' @param param A [cmdstanr::summary()].
#'
#' @returns A [tibble::tibble()].
#' @noRd

clean_summary <- function(summary) {
  id_all <- variable <- NULL
  summary |>
    dplyr::filter(is_learning_model_var(variable)) |>
    dplyr::select(
      variable, mean, tidyselect::any_of(tidyselect::matches("ess|rhat"))
    ) |>
    dplyr::mutate(
      id_all = sub("\\].*$", "", sub(".*\\[", "", variable))
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      id_no = as.numeric(strsplit(id_all, ",")[[1]][1]),
      aff_num = as.numeric(strsplit(id_all, ",")[[1]][2]), # affect models
      outc_lag = as.numeric(strsplit(id_all, ",")[[1]][3]) # affect delta model
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(parameter = sub("\\[.*$", "", variable)) |>
    dplyr::select(-variable, -id_all)
}

#' Construct a tibble containing individual-level parameter values and baseline
#' information
#'
#' \code{make_par_df} combines a [cmdstanr::summary()] with raw data to e.g.,
#' pass to an outcome model, with filtering based on model fit metrics.
#'
#' @param raw Raw data, e.g., saved by [fit_learning_model()].
#' @param summary A [cmdstanr::summary()] data frame.
#' @param rhat_upper Upper bound of split r-hat values to include. Set to
#' \code{Inf} to include all participants.
#' @param ess_lower Lower bound of effective sample size values to include. Set
#' to \code{0} to include all participants.
#' @param bsl_trnsfm A function to apply to intercept. Defaults to identity.
#' @param join_dem Combine output with participant demographic info?
#' @param adj_order Same as [fit_learning_model()].
#'
#' @returns A long format [tibble::tibble()] with model parameters and baseline
#' participant data.
#'
#' @examples \dontrun{
#' fit <- fit_learning_model(
#'   example_data$nd,
#'   model = "2a",
#'   vb = FALSE,
#'   exp_part = "training"
#'  )
#'
#' make_par_df(fit$raw_df, fit$summary, rhat_upper = 1.1, ess_lower = 100)
#' }
#'
#' @export

make_par_df <- function(raw,
                        summary,
                        rhat_upper = Inf,
                        ess_lower = 0,
                        bsl_trnsfm = function(x) x,
                        join_dem = TRUE,
                        adj_order = c("happy", "confident", "engaged")) {

  subjID <- id_no <- aff_num <- parameter <- rhat <- NULL

  if (is.null(raw$subjID)) {
    raw <- raw |> dplyr::mutate(subjID = as.character(id_no))
  }

  ids <- raw |>
    dplyr::distinct(subjID) |>
    dplyr::mutate(id_no = dplyr::row_number())

  n_id <- length(ids$subjID)

  summ <- clean_summary(summary) |>
    dplyr::mutate(
      mean = ifelse(grepl("0", parameter), bsl_trnsfm(mean), mean)
    ) |>
    dplyr::rename(posterior_mean = mean) |>
    dplyr::mutate(adj = ifelse(is.na(aff_num), NA, adj_order[aff_num])) |>
    dplyr::right_join(ids, by = "id_no") |>
    dplyr::group_by(subjID)

  if ("rhat" %in% names(summ) && is.finite(rhat_upper)) {
    summ <- summ |>
      dplyr::filter(!any(rhat > rhat_upper))
  }

  ess_cols <- names(summ)[grepl("^ess_b", names(summ))]
  if (length(ess_cols) > 0 && ess_lower > 0) {
    summ <- summ |>
      dplyr::filter(
        !any(dplyr::if_any(tidyselect::all_of(ess_cols), ~.x < ess_lower))
      )
  }

  summ <- summ |>
    dplyr::ungroup() |>
    dplyr::select(
      tidyselect::vars_select_helpers$where(~!all(is.na(.x)))
    )

  lost_ids <- n_id - length(unique(summ$subjID))
  if (lost_ids > 0) message(
    lost_ids, " individual(s) dropped due to high rhat and/or low bulk ESS."
  )
  if (join_dem) {
    ns_env <- parent.env(environment())

    if (exists("ppt_info", envir = ns_env, inherits = FALSE)) {
      dem <- get("ppt_info", envir = ns_env, inherits = FALSE)
    } else if (is.data.frame(raw)) {
      dem <- raw |>
        dplyr::distinct(subjID, .keep_all = TRUE)
    } else {
      warning(
        "Cannot locate participant demographics; please join manually (e.g.,",
        " from all_res$ppt_info)."
      )
      return(summ)
    }

    summ <- dem |>
      dplyr::inner_join(summ, by = "subjID")
  }
  summ
}

#' Calculate the expected log pointwise predictive density via PSIS-LOO
#'
#' \code{calculate_elpd_loo} computes an approximate leave-one-out
#' cross-validation estimate ([loo::loo()]) from a fitted
#' [cmdstanr::CmdStanMCMC] or [cmdstanr::CmdStanVB] model's \code{log_lik}
#' generated quantity. When \code{fit} comes from variational inference,
#' [loo::loo_approximate_posterior()] is used instead of [loo::loo()], since
#' ordinary PSIS-LOO assumes exact posterior draws.
#'
#' @param fit A [cmdstanr::CmdStanMCMC] or [cmdstanr::CmdStanVB] object, e.g.
#' \code{fit_learning_model(..., outputs = "model_env")$fit}.
#' @param cores Number of cores to use for computing PSIS-LOO. Defaults to
#' \code{options(mc.cores)} or 4 if this is not set.
#'
#' @returns A [loo::loo()] object (class \code{psis_loo} for MCMC fits, or
#' \code{psis_loo_ap} for variational fits).
#'
#' @examples \dontrun{
#' data(example_data)
#' fit <- fit_learning_model(
#'   example_data$nd,
#'   model = "2a",
#'   exp_part = "training",
#'   outputs = "model_env"
#' )
#' calculate_elpd_loo(fit$fit)
#' }
#'
#' @export

calculate_elpd_loo <- function(fit, cores = getOption("mc.cores", 4)) {
  if (!inherits(fit, "CmdStanVB")) {
    return(fit$loo(cores = cores, save_psis = TRUE))
  }

  draws <- fit$draws(
    variables = c("log_lik", "lp__", "lp_approx__"), format = "matrix"
  )
  ll_cols <- grep("^log_lik(\\[|$)", colnames(draws))

  log_lik_mat <- as.matrix(draws[, ll_cols, drop = FALSE])
  log_p <- as.numeric(draws[, "lp__"])
  log_g <- as.numeric(draws[, "lp_approx__"])

  loo::loo_approximate_posterior(
    log_lik_mat, log_p, log_g, cores = cores, save_psis = TRUE
  )
}

#' Check for collapsed between-participant variation in learning-model
#' parameters
#'
#' \code{check_par_sd} flags individual-level parameters (optionally split by
#' affect adjective / outcome lag) whose posterior means show implausibly
#' little variation across participants - a known pathology when fitting
#' complex models (e.g., affect models with many free parameters) with
#' mean-field variational inference. Used internally by
#' [fit_learning_model()] when \code{vb_sd_check = TRUE}, but can also be run
#' directly on any saved fit's summary.
#'
#' @param summary A [cmdstanr::summary()] data frame.
#' @param threshold Minimum acceptable standard deviation of posterior means
#' across participants; parameter/adjective/lag combinations with a lower (or
#' \code{NA}) SD are flagged. Defaults to \code{0.01}.
#' @param ignore_pars Character vector of base parameter names (e.g.
#' \code{"w1_o"}) to exclude from the check, for parameters known to be
#' poorly identified. Defaults to none.
#'
#' @returns A list with \code{sd_df} (standard deviation of posterior means
#' for every parameter/adjective/lag combination, ordered ascending),
#' \code{collapsed} (the subset flagged as below \code{threshold}), and
#' \code{pass} (\code{TRUE} if nothing was flagged).
#'
#' @examples \dontrun{
#' fit <- fit_learning_model(
#'   example_data$nd,
#'   model = "2a",
#'   vb = TRUE,
#'   exp_part = "training",
#'   vb_sd_check = FALSE
#' )
#' check_par_sd(fit$summary)
#' }
#'
#' @export

check_par_sd <- function(summary,
                         threshold = 0.01,
                         ignore_pars = character(0)) {
  parameter <- aff_num <- outc_lag <- mean <- NULL

  sd_df <- clean_summary(summary) |>
    dplyr::filter(!parameter %in% ignore_pars) |>
    dplyr::group_by(parameter, aff_num, outc_lag) |>
    dplyr::summarise(
      sd = stats::sd(mean, na.rm = TRUE), n = dplyr::n(), .groups = "drop"
    ) |>
    dplyr::arrange(sd)

  collapsed <- sd_df |> dplyr::filter(is.na(sd) | sd < threshold)

  list(sd_df = sd_df, collapsed = collapsed, pass = nrow(collapsed) == 0)
}

#' Define axis title name
#'
#' \code{axis_title} returns an axis title, allowing for sub/superscripts and/or
#' Greek letters.
#'
#' @param param Parameter name.
#' @param test,alpha_par Booleans indicating whether the parameter is from a
#' test data and/or is a learning rate parameter.
#' @param alpha_par_nms Names of learning rate parameters (ignored if
#' \code{!alpha_par}).
#' @param mu Boolean indicating if the parameters are group-level.
#'
#' @returns Axis title as a character string containing an expression.
#' @noRd

axis_title <- function(param,
                       p,
                       test,
                       alpha_par,
                       alpha_par_nms,
                       mu = FALSE) {
  if (mu) sub("mu_", "", param)
  if (grepl("w", param, fixed = TRUE)) param <- sub("w", "w_", param)
  spl <- unlist(strsplit(param, "_"))
  s <- ifelse(test, paste0(spl[1], "*minute"), spl[1])
  if (length(spl) != 1) {
    if (alpha_par && !is.null(alpha_par_nms)) {
      a <- paste0("[", alpha_par_nms[p], "]")
    } else {
      a <- paste0("[", spl[2], "]")
    }
    if (length(spl) == 3) {
      sup <- tolower(spl[3])
      if (sup %in% c("neg", "minus")) {
        sup <- '"(-)"'
      } else if (sup %in% c("pos", "plus")) {
        sup <- '"(+)"'
      }
      a <- paste0(a, "^", sup)
    }
  } else {
    a <- ""
  }
  paste0(s, a)
}

#' Extract posterior predictions from affect data models and assess fit
#'
#' \code{get_affect_ppc} combines posterior predictions contained in a
#' [posterior::draws_df()] outputted from a fit model with raw affect ratings,
#' and returns various fit metrics (pseudo-\eqn{R^2}, MAE, RMSE), for each
#' individual. We follow Ferrari & Cribari-Neto (2004) in defining pseudo-R^2 as
#' the squared correlation between observed and mean posterior predictions.
#'
#' @param draws A [posterior::draws_df()]. Draws outputted from
#' [fit_learning_model] (as a [posterior::draws_list()]) should be converted via
#' [posterior::as_draws_df()] - this is memory intensive, hence it is not done
#' internally.
#' @param raw Raw data, e.g., saved by [fit_learning_model()].
#' @param adj Name of the affect adjective - one of "happy", "confident" or
#' "engaged".
#'
#' @returns List containing a dataframe with participant identifiers, numbers,
#' and pseudo \eqn{R^2}, MAE and RMSE for each individual; and a named list (by
#' ID) of data frames with individuals' mean posterior predictions and raw
#' affect data.
#'
#' @examples \dontrun{
#' fit_affect <- fit_learning_model(
#'   example_data$nd,
#'   model = "2a",
#'   affect = TRUE,
#'   exp_part = "training",
#'   algorithm = "fullrank"
#'  )
#'
#'  fit_ls_happy <- get_affect_ppc(
#'    draws = fit_affect$draws,
#'    raw = fit_affect$raw_df,
#'    adj = "happy"
#'  )
#'  }
#'
#' @importFrom data.table .SD as.data.table
#' @export

get_affect_ppc <- function(draws,
                           raw,
                           adj) {

  draws <- as.data.table(draws)
  indiv_ppcs <- list()

  ## to appease R CMD check
  question_type <- trial_no_q <- question_response <- type <- se_pred <-
    "patterns" <- "..aff_tr" <- subjID <- NULL

  n_id <- length(unique(raw$subjID))
  grps <- raw |>
    dplyr::distinct(dplyr::pick(subjID, tidyselect::any_of("distanced")))
  fit_df <- as.data.frame(
    matrix(
      nrow = n_id, ncol = 5,
      dimnames = list(1:n_id, c("subjID", "id_no", "pseudo_R2", "MAE", "RMSE"))
    )
  )
  pb <- txtProgressBar(min = 0, max = n_id, initial = 0, style = 3)

  for (i in 1:n_id) {
    setTxtProgressBar(pb, i)
    id <- unique(raw$subjID)[i]

    affect_raw <- raw |>
      dplyr::filter(question_type == adj & subjID == id) |>
      dplyr::select(trial_no_q, question_response) |>
      dplyr::rename(mean_raw = question_response)

    aff_tr <- which(raw[raw$subjID == id, ]$question_type == adj) #nolint

    affect_pred <-
      draws[, .SD, .SDcols = patterns(paste0("^y_pred\\[", i, ","))]
    affect_pred <- suppressWarnings(affect_pred[, ..aff_tr])
    # https://github.com/Rdatatable/data.table/issues/2988

    means <- colMeans(affect_pred) * 100
    se <- sapply(affect_pred, function(x) sd(x * 100) / sqrt(length(x)))

    pred_df <- tibble::tibble(
      mean_pred = means, se_pred = se
    ) |>
      dplyr::mutate(trial_no_q = dplyr::row_number())

    obs <- affect_raw$mean_raw
    pred <- pred_df$mean_pred

    res_df <- affect_raw |>
      dplyr::left_join(pred_df, by = "trial_no_q") |>
      tidyr::pivot_longer(
        cols = tidyselect::contains("mean"), names_to = "type",
        names_prefix = "mean_"
      ) |>
      dplyr::mutate(se_pred = ifelse(type == "raw", 0, se_pred))

    fit_df[i, 1] <- id
    fit_df[i, 2] <- i
    fit_df[i, 3] <- suppressWarnings(cor(obs, pred)^2) # some have variation = 0
    fit_df[i, 4] <- mean(abs(pred - obs))
    fit_df[i, 5] <- sqrt(mean((pred - obs)^2))

    indiv_ppcs[[i]] <- res_df
  }

  names(indiv_ppcs) <- unique(raw$subjID)

  ret <- list()
  ret$fit_df <- dplyr::left_join(fit_df, grps, by = "subjID")
  ret$indiv_ppcs <- indiv_ppcs

  ret
}
