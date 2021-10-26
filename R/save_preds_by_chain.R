#' Get and save posterior predictions for training data
#'
#' \code{save_preds_by_chain} is a helper function which aims to automate the loading of
#' posterior predictions without crashing R due to memory overload.
#'
#' @param out_files Vector of .csv file names which contain posterior predictions (e.g., outputted from
#' [pstpipeline::generate_posterior_quantities()]).
#' @param out_dir Path to output directory (defaults to current working directory).
#' @param obs_df Raw training data, e.g., outputted from [pstpipeline::fit_learning_model()] (this is best
#' as it ensures individuals are matched with the correct predictions.)
#' @param n_draws_chain Number of MCMC sampling iterations per chain.
#' @param save_outputs Save lists of posterior predictions (per chain) and the updated training data
#' \code{tibble}?
#' @param save_dir Directory to save items to, will be created if it does not exist. Defaults to the
#' directory of the output files.
#' @param prefix Optional prefix to add to the saved objects.
#' @param ... Other arguments which are unlikely to be necessary to change: \code{n_trials} (default = 360);
#' \code{vars} (default = "y_pred"), and
#' \code{pred_types} (default = \code{c("AB", "CD", "EF")}).
#'
#' @return An updated \code{tibble} with summed choices per chain and their overall proportion. This is also
#' saved to disk, plus a named list of posterior predictions for each individual (separately saved by chain).
#'
#' @importFrom magrittr %>%
#' @importFrom rlang := !!
#' @importFrom data.table rbindlist
#' @export

save_preds_by_chain <-
  function(out_files,
           out_dir = "",
           obs_df,
           n_draws_chain,
           save_outputs = TRUE,
           save_dir = out_dir,
           prefix = "",
           ...) {

  start <- Sys.time()

  save_dir <- file.path(getwd(), save_dir)
  if (!dir.exists(save_dir)) dir.create(save_dir)

  n_indiv <- length(unique(obs_df$subjID))
  if (is.null())

  l <- list(...)
  if (is.null(n_trials)) n_trials <- 360
  if (is.null(vars)) vars <- "y_pred"
  if (is.null(pred_types)) pred_types <-  pred_types = c("AB", "CD", "EF")

  paths <- out_files
  if (out_dir != "") {
    for (o in seq_along(out_files)) {
      paths[o] <- paste0(out_dir, out_files[o])
    }
  }

  var_order <- as.vector(
    sapply(1:n_indiv, function(i)
    sapply(1:n_trials, function(j) paste0("y_pred[", i, ",", j, "]"))
    )
  )

  indiv <- unique(obs_df$subjID)
  indiv_obs_df <- data.table::data.table()

  for (o in seq_along(paths)) {
    indiv_obs <- list()

    indiv_draws <- cmdstanr::read_cmdstan_csv(paths[o], variables = vars, format = "draws_list")[[2]][[1]] %>%
      dplyr::select(tidyselect::all_of(var_order)) %>%
      split.default(rep(1:n_indiv, each = n_trials))
    indiv_draws <- lapply(names(indiv_draws), FUN = function(d) indiv_draws[[d]][, colSums(indiv_draws[[d]]) >= 0])
      # as columns missing PPCs are set to -1

    for (i in seq_along(indiv_draws)) {
      pred_indiv <- indiv_draws[[i]]
      missing <- which(!seq(1:n_trials) %in% obs_df[obs_df$subjID == indiv[i],]$trial_no)
      for (m in missing) {
        pred_indiv <- pred_indiv %>%
          tibble::add_column(NA, .after = m-1)
      }
      colnames(pred_indiv) <- sapply(1:n_trials, function(n) paste0("y_pred_", n))

      indiv_obs_types <- list()

      nwnme = rlang::sym(paste0("choice_pred_sum_ch_", o))

      if (any(pred_types == "AB")) {
        ## add observed/predicted choices at correct indexes
        ab_pred_nms <- paste0("y_pred_", obs_df[obs_df$subjID == indiv[i] & obs_df$type == 12,]$trial_no)
        ab_preds <- pred_indiv %>%
          dplyr::select(tidyselect::all_of(ab_pred_nms))

        ab_pred_sums <- tibble::as_tibble(colSums(ab_preds), rownames = "trial_no") %>%
          dplyr::rename(!!nwnme := value) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(trial_no = as.integer(tail(strsplit(trial_no, "_")[[1]], n=1))) %>%
          dplyr::ungroup()

        ab_obs <- tibble::as_tibble(obs_df[obs_df$subjID == indiv[i] & obs_df$type == 12]) %>%
          dplyr::mutate(id_no = i) %>%
          dplyr::left_join(ab_pred_sums, by = "trial_no")

        indiv_obs_types$ab <- ab_obs
      }

      if (any(pred_types == "CD")) {
        ## add observed/predicted choices at correct indexes
        cd_pred_nms <- paste0("y_pred_", obs_df[obs_df$subjID == indiv[i] & obs_df$type == 34,]$trial_no)
        cd_preds <- pred_indiv %>%
          dplyr::select(tidyselect::all_of(cd_pred_nms))

        cd_pred_sums <- tibble::as_tibble(colSums(cd_preds), rownames = "trial_no") %>%
          dplyr::rename(!!nwnme := value) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(trial_no = as.integer(tail(strsplit(trial_no, "_")[[1]], n=1))) %>%
          dplyr::ungroup()

        cd_obs <- tibble::as_tibble(obs_df[obs_df$subjID == indiv[i] & obs_df$type == 34]) %>%
          dplyr::mutate(id_no = i) %>%
          dplyr::left_join(cd_pred_sums, by = "trial_no")

        indiv_obs_types$cd <- cd_obs
      }
      if (any(pred_types == "EF")) {
        ef_pred_nms <- paste0("y_pred_", obs_df[obs_df$subjID == indiv[i] & obs_df$type == 56,]$trial_no)
        ef_preds <- pred_indiv %>%
          dplyr::select(tidyselect::all_of(ef_pred_nms))

        ef_pred_sums <- tibble::as_tibble(colSums(ef_preds), rownames = "trial_no") %>%
          dplyr::rename(!!nwnme := value) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(trial_no = as.integer(tail(strsplit(trial_no, "_")[[1]], n=1))) %>%
          dplyr::ungroup()

        ef_obs <- tibble::as_tibble(obs_df[obs_df$subjID == indiv[i] & obs_df$type == 56]) %>%
          dplyr::mutate(id_no = i) %>%
          dplyr::left_join(ef_pred_sums, by = "trial_no")

        indiv_obs_types$ef <- ef_obs
      }
      indiv_draws[[i]] <- pred_indiv
      indiv_obs[[i]] <- data.table::rbindlist(indiv_obs_types)
    }

    names(indiv_draws) <- indiv
    if (save_outputs) saveRDS(indiv_draws, paste0(save_dir, "/", prefix, "chain_", o, "_ppc_list.RDS"))
    rm(indiv_draws)

    if (length(indiv_obs_df) == 0) {
      indiv_obs_df <- data.table::rbindlist(indiv_obs)
    } else {
      indiv_df_oc <- data.table::rbindlist(indiv_obs) %>%
        dplyr::select(subjID, trial_no, paste0("choice_pred_sum_ch_", o))
      indiv_obs_df <- indiv_obs_df %>%
        dplyr::left_join(indiv_df_oc, by = c("subjID", "trial_no"))
    }

    if (o == length(paths)) {
      indiv_obs_df <- indiv_obs_df %>%
        dplyr::rowwise() %>%
        dplyr::mutate(choice_pred_prop = sum(
          dplyr::across(tidyselect::contains("pred_sum")))/(n_draws_chain*o)
          ) %>%
        dplyr::ungroup()
    }

    rm(indiv_obs)
  }

  if (save_outputs) saveRDS(indiv_obs_df, file = paste0(save_dir, "/", prefix, "indiv_obs_sum_ppcs_df.RDS"))
  message(paste0("Finished in ", round(difftime(Sys.time(), start, unit = "secs")[[1]], digits = 1), " seconds."))

  return(indiv_obs_df)
}
