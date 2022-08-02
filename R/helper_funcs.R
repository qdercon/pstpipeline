#' @noRd
#' @keywords internal

single_hdi <- function(vars,
                       cred) {
  # adapted from hBayesDM::HDIofMCMC
  sampleVec <- as.vector(t(vars))
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
      HDI_lower <- single_hdi(vars = var, cred = cred_mass)[1]
      returns <- cbind(returns, HDI_lower)
    } else {
      cred_mass = 2*q - 1
      HDI_upper <- single_hdi(vars = var, cred = cred_mass)[2]
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
family_ch <- function(param) {
  if (grepl("alpha|gamma", param)) return(Gamma(link = "log"))
  else return(gaussian())
}
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
axis_title <- function(param, p, test, alpha_par, alpha_par_nms) {
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
