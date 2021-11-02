#' Get and store posterior predictions for training data
#'
#' \code{get_preds_by_chain} is a helper function which aims to automate the loading of posterior predictions
#' without crashing R due to memory overload. This is done simply by importing chain-by-chain predictions, and
#' summing over all chains/draws the predicted binary choices for each individual and trial.
#'
#' @param out_files Vector of .csv file names which contain posterior predictions (e.g., outputted from
#' [pstpipeline::generate_posterior_quantities()]).
#' @param out_dir Path to output directory (defaults to current working directory).
#' @param obs_df Raw training data, e.g., outputted from [pstpipeline::fit_learning_model()] (this is best
#' as it ensures individuals are matched with the correct predictions.)
#' @param n_draws_chain Number of MCMC sampling iterations per chain.
#' @param save_dir Directory to save items to, will be created if it does not exist. Defaults to the
#' directory of the output files.
#' @param prefix Optional prefix to add to the saved objects.
#' @param memory_save An alternative method to obtain predictions, which loads the predictions for each
#' individual (across all chains) one-by-one, as opposed to importing all the draws for all individuals.
#' This will be significantly slower but enables the function to run with very limited RAM.
#' @param ... Other arguments which are unlikely to be necessary to change: \code{n_trials} (default = 360);
#' \code{vars} (default = "y_pred"), and \code{pred_types} (default = \code{c("AB", "CD", "EF")}).
#'
#' @return An updated \code{tibble} with summed choices per chain and their overall proportion.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang := !!
#' @importFrom data.table rbindlist
#' @export

get_preds_by_chain <-
  function(out_files,
           out_dir = "",
           obs_df,
           n_draws_chain,
           save_dir = out_dir,
           test = FALSE,
           prefix = "",
           splits = list(blocks = 1:6, sum_blks = list(c(1,3), c(4,6))),
           exclude = NULL,
           memory_save = TRUE,
           ...) {

  start <- Sys.time()

  save_dir <- file.path(getwd(), save_dir)
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  l <- list(...)
  if (!test) {
    if (is.null(l$n_trials)) l$n_trials <- 360
    if (is.null(l$n_blocks)) l$n_blocks <- 6
    if (is.null(l$pred_var)) l$pred_var <- "y_pred"
    if (is.null(l$pred_types)) l$pred_types <- c("AB", "CD", "EF")
    n_indiv <- length(unique(obs_df$subjID))
  }
  else {
    all_pairs <- list("12", "34", "56", "13", "14", "15", "16", "32", "42",
                      "52", "62","35", "36", "54", "64")
    names(all_pairs) <- c("AB", "CD", "EF", "AC", "AD", "AE", "AF", "CB",
                          "DB", "EB", "FB", "CE", "CF", "ED", "FD")

    if (is.null(l$n_trials)) l$n_trials <- 60
    if (is.null(l$n_blocks)) l$n_blocks <- 1
    if (is.null(l$pred_var)) l$pred_var <- "y_pred"
    if (is.null(l$pred_types)) l$pred_types <- names(all_pairs)
    splits <- list(list(), list())
    obs_df <- obs_df$test
    n_indiv <- length(unique(obs_df$subjID))
  }

  trials_per_block <- l$n_trials/l$n_blocks

  paths <- out_files
  if (out_dir != "") {
    for (o in seq_along(out_files)) {
      paths[o] <- paste0(out_dir, "/", out_files[o])
    }
  }

  pred_var_names <- as.vector(
    sapply(1:n_indiv, function(i)
    sapply(1:l$n_trials, function(j) paste0(l$pred_var, "[", i, ",", j, "]"))
    )
  )
  indiv_vars <- split(pred_var_names, ceiling(seq_along(pred_var_names)/l$n_trials))

  indiv <- unique(obs_df$subjID)
  ids <- sapply(1:n_indiv, function(i) list(i))
  names(ids) <- indiv

  indiv_obs_list <- obs_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(id_no = ids[[subjID]]) %>%
    dplyr::ungroup() %>%
    split(., as.factor(.$id_no))

  trial_avg_types <- 1+sum(length(splits[[1]]), length(splits[[2]]))
  trial_avg_df <- tibble::as_tibble(
    matrix(
      ncol = 8,
      nrow = 0,
      dimnames = list(NULL, c("subjID", "id_no", "obs_mean", "pred_post_mean",
                              "pred_post_median", "pred_post_lower_95_hdi",
                              "pred_post_upper_95_hdi", "type"))
    )
  ) %>%
    dplyr::mutate(dplyr::across(c(1,8), as.character)) %>%
    dplyr::mutate(dplyr::across(2:7, as.numeric))
  trial_avg_list <- rep(list(trial_avg_df), n_indiv)

  if (!memory_save) {
    all_draws <- list()
    for (o in seq_along(paths)) {
      all_draws[[o]] <-
        cmdstanr::read_cmdstan_csv(paths[o], variables = l$pred_var, format = "draws_list")[[2]][[1]]
      if (is.list(all_draws[[o]])) all_draws[[o]] <- data.table::as.data.table(all_draws[[o]])
    }
    all_draws_df <- data.table::rbindlist(all_draws)
    rm(all_draws)
  }

  pb = txtProgressBar(min = 0, max = n_indiv, initial = 0, style = 3)

  # first get the trial-wise summed predictions, by individual (i.e., the proportion of choices for trial 1, 2, etc...)
  # taking care to line up predictions correctly - these will then be aggregated across individuals

  for (id in seq_along(indiv)) {
    setTxtProgressBar(pb, id)
    if (memory_save) {
      all_indiv_draws <- list()
      for (o in seq_along(paths)) {
        all_indiv_draws[[o]] <-
          cmdstanr::read_cmdstan_csv(paths[o], variables = indiv_vars[[id]], format = "draws_list")[[2]][[1]]
        if (is.list(all_indiv_draws[[o]])) all_indiv_draws[[o]] <- data.table::as.data.table(all_indiv_draws[[o]])
        all_indiv_draws[[o]] <- all_indiv_draws[[o]] %>% dplyr::select(-where(~ any(. == -1)))
      }
      all_indiv_draws <- data.table::rbindlist(all_indiv_draws)
    }
    else {
      all_indiv_draws <- all_draws_df %>%
        dplyr::select(indiv_vars[[id]]) %>%
        dplyr::select(-where(~ any(. == -1)))
    }

    indiv_obs_list_id <- indiv_obs_list[[id]]

    missing <- which(!seq(1:l$n_trials) %in% indiv_obs_list_id$trial_no)
    for (m in missing) {
      all_indiv_draws <- all_indiv_draws %>%
        tibble::add_column(NA, .after = m-1) # adds columns of NAs in the correct place where a trial was missed
    }
    colnames(all_indiv_draws) <- sapply(1:l$n_trials, function(n) paste(l$pred_var, n, sep = "_"))

    indiv_obs_types <- list()

    if (test) type_key <- all_pairs
    else {
      type_key <- list(12, 34, 56)
      names(type_key) <- c("AB", "CD", "EF")
    }
    for (stim in seq_along(l$pred_types)) {
      ## add observed/predicted choices at correct indexes
      pred_nms <- paste0(
        l$pred_var, "_", indiv_obs_list_id[indiv_obs_list_id$type == type_key[[l$pred_types[stim]]],]$trial_no
        )
      preds <- all_indiv_draws %>%
        dplyr::select(tidyselect::all_of(pred_nms))

      pred_sums <- tibble::as_tibble(colSums(preds)/dim(preds)[1], rownames = "trial_no") %>%
        dplyr::rename(choice_pred_prop = value) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(trial_no = as.integer(tail(strsplit(trial_no, "_")[[1]], n=1))) %>%
        dplyr::ungroup()

      obs <-
        tibble::as_tibble(
          indiv_obs_list_id[indiv_obs_list_id$type == type_key[[l$pred_types[stim]]],]
        ) %>%
        dplyr::mutate(id_no = id) %>%
        dplyr::left_join(pred_sums, by = "trial_no")

      indiv_obs_types[[l$pred_types[stim]]] <- obs

      ## get different types of rowsum
      summed_trials <- rowSums(preds)/dim(preds)[2]
      all_trials_q <- pstpipeline::quantile_hdi(summed_trials, c(0.025, 0.5, 0.975))
      trial_avg_list[[id]] <- trial_avg_list[[id]] %>%
        dplyr::bind_rows(
          tibble::tibble(
            "subjID" = indiv[id],
            "id_no" = id,
            "obs_mean" = mean(obs$choice),
            "pred_post_mean" = mean(summed_trials),
            "pred_post_median" = all_trials_q[2],
            "pred_post_lower_95_hdi" = all_trials_q[1],
            "pred_post_upper_95_hdi" = all_trials_q[3],
            "type" = paste(l$pred_types[stim], "all_trials", sep = "_")
          )
        )

      if (length(splits[[1]]) > 0 | length(splits[[2]]) > 0) {
        max_trial <- seq(0, l$n_trials, trials_per_block)
        trial_nums <- sapply(strsplit(names(preds), "_"), function(g) return(as.integer(g[3])))
        pred_names <- lapply(
          split(trial_nums, cut(trial_nums, max_trial)), function(h) h <- paste0(l$pred_var, "_", h)
        )
        names(pred_names) <- paste0("block_", 1:l$n_blocks)
        blknames <- vector("character")

        blk_indiv <- length(splits[[1]])
        grpd <- length(splits[[2]])

        ## individual blocks
        if (length(blk_indiv) > 0) {
          blknames <- paste0("block_", splits[[1]])
          for (blk in seq_along(blknames)) {
            preds_blk <- preds %>%
              dplyr::select(tidyselect::all_of(pred_names[[blknames[blk]]]))
            summed_trials_blk <- rowSums(preds_blk)/dim(preds_blk)[2]
            blk_trials_q <- pstpipeline::quantile_hdi(summed_trials_blk, c(0.025, 0.5, 0.975))

            trial_avg_list[[id]] <- trial_avg_list[[id]] %>%
              dplyr::bind_rows(
                tibble::tibble(
                  "subjID" = indiv[id],
                  "id_no" = id,
                  "obs_mean" = mean(obs[obs$trial_no > ((blk-1)*trials_per_block) &
                                        obs$trial_no <= (blk*trials_per_block),]$choice),
                  "pred_post_mean" = mean(summed_trials_blk),
                  "pred_post_median" = blk_trials_q[2],
                  "pred_post_lower_95_hdi" = blk_trials_q[1],
                  "pred_post_upper_95_hdi" = blk_trials_q[3],
                  "type" = paste(l$pred_types[stim], blknames[blk], sep = "_")
                  )
              )
          }
        }
        ## grouped blocks
        if (grpd > 0) {
          blkgrpnms <- lapply(
            1:length(splits[[2]]),
            function(b) paste("block", splits[[2]][[b]][1], "to", splits[[2]][[b]][2], sep = "_")
          )

          for (blkgrp in seq_along(blkgrpnms)) {
            included <- seq(splits[[2]][[blkgrp]][1], splits[[2]][[blkgrp]][2])
            included_names <- as.vector(unlist(sapply(included, function(nm) pred_names[[blknames[nm]]])))
            preds_blkgrp <- preds %>%
              dplyr::select(tidyselect::all_of(included_names))
            summed_trials_blkgrp <- rowSums(preds_blkgrp)/dim(preds_blkgrp)[2]
            blk_trial_grp_q <- pstpipeline::quantile_hdi(summed_trials_blkgrp, c(0.025, 0.5, 0.975))

            trial_avg_list[[id]] <- trial_avg_list[[id]] %>%
              dplyr::bind_rows(
                tibble::tibble(
                  "subjID" = indiv[id],
                  "id_no" = id,
                  "obs_mean" = mean(obs[obs$trial_no > (min(included-1)*trials_per_block) &
                                    obs$trial_no <= (max(included)*trials_per_block),]$choice),
                  "pred_post_mean" = mean(summed_trials_blkgrp),
                  "pred_post_median" = blk_trial_grp_q[2],
                  "pred_post_lower_95_hdi" = blk_trial_grp_q[1],
                  "pred_post_upper_95_hdi" = blk_trial_grp_q[3],
                  "type" = paste(l$pred_types[stim], blkgrpnms[blkgrp], sep = "_")
                )
              )
          }
        }
      }
    }
    if (test) {
      test_grps <- c("training", "chooseA", "avoidB", "novel")
      for (test_grp in test_grps) {
        pred_nms <- paste0(
          l$pred_var, "_", indiv_obs_list_id[indiv_obs_list_id$test_type == test_grp,]$trial_no
        )
        preds <- all_indiv_draws %>%
          dplyr::select(tidyselect::all_of(pred_nms))
        obs <-
          tibble::as_tibble(
            indiv_obs_list_id[indiv_obs_list_id$type == type_key[[l$pred_types[stim]]],]
          )
        group_sum_trials <- rowSums(preds)/dim(preds)[2]
        group_trials_q <- pstpipeline::quantile_hdi(group_sum_trials, c(0.025, 0.5, 0.975))
        trial_avg_list[[id]] <- trial_avg_list[[id]] %>%
          dplyr::bind_rows(
            tibble::tibble(
              "subjID" = indiv[id],
              "id_no" = id,
              "obs_mean" = mean(obs$choice),
              "pred_post_mean" = mean(group_sum_trials),
              "pred_post_median" = group_trials_q[2],
              "pred_post_lower_95_hdi" = group_trials_q[1],
              "pred_post_upper_95_hdi" = group_trials_q[3],
              "type" = paste(test_grp, "all_trials", sep = "_")
            )
          )
      }
    }
    indiv_obs_list[[id]] <- data.table::rbindlist(indiv_obs_types)
  }

  indiv_obs_df <- data.table::rbindlist(indiv_obs_list)
  trial_obs_df <- data.table::rbindlist(trial_avg_list)

  if (!is.null(exclude)) {
    indiv_obs_df <- indiv_obs_df %>% dplyr::filter(!id_no %in% exclude)
    trial_obs_df <- trial_obs_df %>% dplyr::filter(!id_no %in% exclude)
  }

  saveRDS(indiv_obs_df, file = paste0(save_dir, "/", prefix, "indiv_obs_sum_ppcs_df.RDS"))
  saveRDS(trial_obs_df, file = paste0(save_dir, "/", prefix, "trial_block_avg_hdi_ppcs_df.RDS"))

  message(paste0("Finished in ", round(difftime(Sys.time(), start, unit = "secs")[[1]], digits = 1), " seconds."))

  ret <- list()
  ret$indiv_obs_df <- indiv_obs_df
  ret$trial_obs_df <- trial_obs_df

  return(ret)
}
