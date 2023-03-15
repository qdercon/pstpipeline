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

  return(subsample)
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
  return(as.vector(t(HDIlim)))
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
                         ...) {

  l <- list(...)
  if (is.null(l$transform)) l$transform <- FALSE
  if (l$transform) { # fixes issues calling function from within ggplot with
                     # exponentiated coefficients
    var <- log(var / 100 + 1)
  }
  returns <- vector(mode = "numeric")

  for (q in quantile) {
    if (q == 0.5) {
      returns <- cbind(returns, median(var))
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

  return(ret)
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
  if (grepl("alpha|gamma", param)) return(Gamma(link = "log"))
  else return(gaussian())
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
    dplyr::filter(grepl("alpha|beta|w|gamma\\[", variable)) |>
    dplyr::filter(!grepl("_pr|_s|mu|sigma|_diff", variable)) |>
    dplyr::select(
      variable, mean, tidyselect::any_of(tidyselect::matches("ess|rhat"))
    ) |>
    dplyr::mutate(
      id_all = sub("\\].*$", "", sub(".*\\[", "", variable))
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      id_no = as.numeric(strsplit(id_all, ",")[[1]][1]),
      aff_num = as.numeric(strsplit(id_all, ",")[[1]][2])
    ) |>
    # NA unless affect model
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
                        rhat_upper,
                        ess_lower,
                        join_dem = TRUE,
                        adj_order = c("happy", "confident", "engaged")) {

  subjID <- aff_num <- NULL

  ids <- raw |>
    dplyr::distinct(subjID) |>
    dplyr::mutate(id_no = dplyr::row_number())

  n_id <- length(ids$subjID)

  summ <- clean_summary(summary) |>
    dplyr::rename(posterior_mean = mean) |>
    dplyr::mutate(adj = ifelse(is.na(aff_num), NA, adj_order[aff_num])) |>
    dplyr::right_join(ids, by = "id_no") |>
    dplyr::group_by(subjID) |>
    dplyr::filter(dplyr::if_any(
      tidyselect::any_of("rhat"), ~!any(.x > rhat_upper))
    ) |>
    dplyr::filter(dplyr::if_any(
      tidyselect::any_of(
        tidyselect::starts_with("ess_b")), ~!any(.x < ess_lower)
    )) |>
    dplyr::select(
      tidyselect::vars_select_helpers$where(~!all(is.na(.x)))
    ) |>
    dplyr::ungroup()

  lost_ids <- n_id - length(unique(summ$subjID))
  if (lost_ids > 0) message(
    lost_ids, " individual(s) dropped due to high rhat and/or low bulk ESS."
    )
  if (join_dem) {
    summ <- ppt_info |>
      dplyr::inner_join(summ, by = "subjID")
  }
  return(summ)
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
      a <- paste0(a, "^", spl[3])
    }
  } else {
    a <- ""
  }
  return(paste0(s, a))
}

#' Extract posterior predictions from affect data models and assess fit
#'
#' \code{get_affect_ppc} combines posterior predictions contained in a
#' [posterior::draws_df()] outputted from a fit model with raw affect ratings,
#' and returns various fit metrics (\eqn{R^2}, MAE, RMSE), for each individual.
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
#' and \eqn{R^2}, MAE and RMSE for each individual; and a named list (by ID) of
#' data frames with individuals' mean posterior predictions and raw affect data.
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
    matrix(nrow = n_id, ncol = 5,
           dimnames = list(1:n_id, c("subjID", "id_no", "R2", "MAE", "RMSE")))
  )
  pb <- txtProgressBar(min = 0, max = n_id, initial = 0, style = 3)

  for (i in 1:n_id) {
    setTxtProgressBar(pb, i)
    id <- unique(raw$subjID)[i]

    affect_raw <- raw |>
      dplyr::filter(question_type == adj & subjID == id) |>
      dplyr::select(trial_no_q, question_response) |>
      dplyr::rename(mean_raw = question_response)

    aff_tr <- which(raw[raw$subjID == id,]$question_type == adj) #nolint

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

  return(ret)
}
