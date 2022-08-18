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
#' @return A [base::vector()] with the lower and upper HDI.
#' @noRd

# Adapted from hBayesDM::HDIofMCMC, in turn based on John Kruschke's code.

single_hdi <- function(vals,
                       cred) {
  # adapted from
  sampleVec <- as.vector(t(vals))
  sortedPts = sort(sampleVec)
  ciIdxInc = floor(cred * length(sortedPts))
  nCIs = length(sortedPts) - ciIdxInc
  ciWidth = rep(0 , nCIs)
  for (i in 1:nCIs) {
    ciWidth[i] = sortedPts[i + ciIdxInc] - sortedPts[i]
  }
  HDImin = sortedPts[which.min(ciWidth)]
  HDImax = sortedPts[which.min(ciWidth) + ciIdxInc]
  HDIlim = c(HDImin , HDImax)
  return(as.vector(t(HDIlim)))
}

#' Compute quantiles of a probability distrbution based on highest density
#' intervals (HDIs)
#'
#' \code{quantile_hdi} computes any number of highest density quantiles from a
#' sample of representative values, estimated as shortest credible intervals.
#'
#' @param var A vector of representative values from a probability distribution
#' (e.g., MCMC samples).
#' @param quantile A vector of quantiles to return.
#' @param transform Are values log-scaled (e.g., estimated from Gamma GLMs)?
#'
#' @return A sorted [base::vector()] with all HDIs plus the min/max/median if
#' specified.
#' @export

quantile_hdi <- function(var,
                         quantile,
                         transform = FALSE) {

  if (transform) {
   var <- log(var/100 + 1)
  }

  returns <- vector(mode = "numeric")
  for (q in quantile) {
    if (q == 0.5) {
      returns <- cbind(returns, median(var))
    } else if (q == 0) {
      returns <- cbind(returns, min(var))
    } else if (q == 1){
      returns <- cbind(returns, max(var))
    } else if (q < 0.5) {
      cred_mass = 1 - 2*q
      HDI_lower <- single_hdi(vals = var, cred = cred_mass)[1]
      returns <- cbind(returns, HDI_lower)
    } else {
      cred_mass = 2*q - 1
      HDI_upper <- single_hdi(vals = var, cred = cred_mass)[2]
      returns <- cbind(returns, HDI_upper)
    }
  }

  ret <- sort(returns)
  if (transform) {
    ret <- (exp(ret) - 1) *100
  }
  names(ret) <- sapply(1:length(quantile),
                       FUN = function (x) paste0(quantile[x] * 100, "%"))

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
#' @return A [stats::family()].
#' @noRd

family_ch <- function(param) {
  if (grepl("alpha|gamma", param)) return(Gamma(link = "log"))
  else return(gaussian())
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
#'
#' @return A long format [tibble::tibble()] with model parameters and baseline
#' participant data.
#' @importFrom magrittr %>%
#' @export

make_par_df <- function(raw,
                        summary,
                        rhat_upper,
                        ess_lower) {
  subjID <- variable <- . <- matches <- NULL # appease R CMD check
  ids <- raw %>%
    dplyr::distinct(subjID, .keep_all = TRUE) %>%
    dplyr::mutate(id_no = dplyr::row_number())

  n_id <- length(ids$subjID)

  summ <- summary %>%
    dplyr::filter(grepl("alpha|beta|gamma|w", variable)) %>%
    dplyr::filter(!grepl("_pr|mu_|sigma|_s", variable)) %>%
    dplyr::select(
      variable, mean, tidyselect::any_of(matches("ess|rhat"))
    ) %>%
    dplyr::mutate(
      id_all = sub("\\].*$", "", sub(".*\\[", "", .[["variable"]]))
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      id_no = as.numeric(strsplit(id_all, ",")[[1]][1]),
      aff_num = as.numeric(strsplit(id_all, ",")[[1]][2])
    ) %>%
    # NA unless affect model
    dplyr::ungroup() %>%
    dplyr::mutate(parameter = sub("\\[.*$", "", .[["variable"]])) %>%
    dplyr::select(-variable, -id_all) %>%
    dplyr::rename(posterior_mean = mean) %>%
    dplyr::right_join(ids, by = "id_no") %>%
    dplyr::group_by(subjID) %>%
    dplyr::filter(dplyr::if_any(
      tidyselect::any_of("rhat"), ~!any(.x > rhat_upper))
    ) %>%
    dplyr::filter(dplyr::if_any(
      tidyselect::any_of(
        tidyselect::starts_with("ess_b")), ~!any(.x < ess_lower)
    )) %>%
    dplyr::select(
      tidyselect::vars_select_helpers$where(~!all(is.na(.x)))
    ) %>%
    dplyr::ungroup()

  lost_ids <- n_id - length(unique(summ$subjID))
  if (lost_ids > 0) message(
    lost_ids, " individual(s) dropped due to high rhat and/or low bulk ESS."
    )

  par_df <- ppt_info %>%
    dplyr::inner_join(summ, by = "subjID")

  return(par_df)
}

#' Define axis title name
#'
#' \code{axis_title} returns an axis title, allowing for sub/superscripts and/or
#' Greek letters.
#'
#' @param param Parameter name.
#' @param p Iteration number.
#' @param test,alpha_par Booleans indicating whether the parameter is from a
#' test data and/or is a learning rate parameter.
#' @param alpha_par_nms Names of learning rate parameters (ignored if
#' \code{!alpha_par}).
#'
#' @return Axis title as a character string containing an expression.
#' @noRd

axis_title <- function(param,
                       p,
                       test,
                       alpha_par,
                       alpha_par_nms) {
  if (grepl("w", param, fixed = TRUE)) param <- sub("w", "w_", param)
  spl <- unlist(strsplit(param, "_"))
  s <- ifelse(test, paste0(spl[1], "*minute"), spl[1])
  if (length(spl) != 1) {
    if (alpha_par & !is.null(alpha_par_nms)) {
      a <- paste0("[", alpha_par_nms[p], "]")
    }
    else a <- paste0("[", spl[2], "]")
    if (length(spl) == 3) {
      a <- paste0(a, "^", spl[3])
    }
  } else a <- ""
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
#' @return List containing a [base::data.frame()] with participant identifiers,
#' numbers, and \eqn{R^2}, MAE and RMSE for each individual; and a named list
#' (by ID) of data frames with individuals' mean posterior predictions and raw
#' affect data.
#' @importFrom magrittr %>%
#' @importFrom data.table .SD as.data.table
#' @export

get_affect_ppc <- function(draws,
                           raw,
                           adj) {

  draws <- as.data.table(draws)
  indiv_ppcs <- list()

  n_id <- length(unique(raw$subjID))
  fit_df <- as.data.frame(
    matrix(nrow = n_id, ncol = 5,
           dimnames = list(1:n_id, c("subjID", "id_no", "R2", "MAE", "RMSE")))
  )
  pb = txtProgressBar(min = 0, max = n_id, initial = 0, style = 3)

  for (i in 1:n_id) {
    setTxtProgressBar(pb, i)
    id <- unique(raw$subjID)[i]

    affect_raw <- raw %>%
      dplyr::filter(question_type == adj & subjID == id) %>%
      dplyr::select(trial_no_q, question_response) %>%
      dplyr::rename(mean_raw = question_response)

    aff_tr <- which(raw[raw$subjID == id,]$question_type == adj)

    affect_pred <-
      draws[, .SD, .SDcols = patterns(paste0("^y_pred\\[", i, ","))]
    affect_pred <- affect_pred[, ..aff_tr]

    means <- colMeans(affect_pred)*100
    se <- sapply(affect_pred, function(x) sd(x*100)/sqrt(length(x)))

    pred_df <- tibble::tibble(
      mean_pred = means, se_pred = se
    ) %>%
      dplyr::mutate(trial_no_q = dplyr::row_number())

    obs <- affect_raw$mean_raw
    pred <- pred_df$mean_pred

    res_df <- affect_raw %>%
      dplyr::left_join(pred_df, by = "trial_no_q") %>%
      tidyr::pivot_longer(
        cols = tidyselect::contains("mean"), names_to = "type",
        names_prefix = "mean_"
      )

    fit_df[i,1] <- id
    fit_df[i,2] <- i
    fit_df[i,3] <- suppressWarnings(cor(obs, pred)^2) # some have variation = 0
    fit_df[i,4] <- mean(abs(pred-obs))
    fit_df[i,5] <- sqrt(mean((pred-obs)^2))

    indiv_ppcs[[i]] <- res_df
  }

  names(indiv_ppcs) <- unique(raw$subjID)

  ret <- list()
  ret$fit_df <- fit_df
  ret$indiv_ppcs <- indiv_ppcs

  return(ret)
}

#' Extract individual-level weight parameters from affect data models
#'
#' \code{get_affect_wts} obtains weight parameters from affect data models, and
#' names them with the appropriate affect adjective.
#'
#' @param summary A [cmdstanr::summary()].
#' @param adj_order Same as [fit_learning_model()].
#'
#' @return A [tibble::tibble()] with weights by individual and adjective.
#' @importFrom magrittr %>%
#' @export

get_affect_wts <- function(summary,
                           adj_order = c("happy", "confident", "engaged")) {

  summary %>%
    dplyr::filter(grepl("w|gamma\\[", variable)) %>%
    dplyr::filter(!grepl("mu|sigma", variable)) %>%
    dplyr::select(variable, mean) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      par = unlist(strsplit(gsub("\\].*", "", variable), "\\[|,"))[1],
      id_no = as.numeric(unlist(strsplit(gsub("\\].*", "", variable), "\\[|,"))[2]),
      aff_num = as.numeric(unlist(strsplit(gsub("\\].*", "", variable), "\\[|,"))[3]),
      adj = adj_order[aff_num]
    ) %>%
    dplyr::rename(post_mean = mean) %>%
    dplyr::ungroup() %>%
    dplyr::select(id_no, par, aff_num, post_mean, adj)
}
