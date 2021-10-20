#' Load posterior predictions for training data
#'
#' \code{load_posterior_quantities} is a helper function which aims to automate the loading of
#' posterior predictions without crashing R due to memory overload.
#'
#' @param out_files csvs
#' @param out_dir path
#'
#' @return Series of .RDS files with averaged predictions at different levels (saved to disk).
#'
#' @importFrom magrittr %>%
#' @export

load_posterior_quantities <-
  function(out_files,
           obs_df,
           n_draws_chain,
           out_dir = "",
           n_indiv = NULL,
           prefix = "",
           n_trials = 360,
           sub_dir = "",
           vars = "y_pred",
           split_indiv = FALSE,
           obs = NULL,
           pred_type = c("AB", "CD", "EF"),
           n_chains = NULL,
           chains_saved = FALSE,
           ...) {

  start <- Sys.time()

  save_dir <- file.path(getwd(), sub_dir)
  if (!dir.exists(save_dir)) dir.create(save_dir)

  if (is.null(n_indiv)) n_indiv <- length(unique(obs_df$subjID))

  l <- list(...)

  if (is.null(l$reorder)) l$reorder <- TRUE
  if (is.null(l$separated_indiv)) l$separated_indiv <- split_indiv
  if (is.null(l$ret_obs)) l$ret_obs <- TRUE
    # no need to add confusion!

  paths <- out_files
  if (out_dir != "") {
    for (o in seq_along(out_files)) {
      paths[o] <- paste0(out_dir, out_files[o])
    }
  }

  if (!chains_saved) {
    import_save_chains(
      out_files = paths, n_indiv = n_indiv, n_trials = n_trials, prefix = prefix, reorder = l$reorder,
      split_indiv = split_indiv, sub_dir = save_dir
    )
  }

  for (ch in 1:n_chains) {
    save_obs <- ifelse(ch == 1 & l$ret_obs, TRUE, FALSE)

    plot_list <- pstpipeline::prep_data(
      obs = obs_df, pred_type = pred_type, n_trials = n_trials, n_draws = n_draws, prefix = prefix,
      separated_indiv = split_indiv, sub_dir = save_dir, chain = ch, ret_obs = save_obs
    )

    if (save_obs) saveRDS(plot_list$obs_list, file = paste0(save_dir, "/observed.RDS"))
    for (p in pred_type) {
      saveRDS(plot_list$pred_list[[p]],
              file = paste0(save_dir, "/chain", ch, "_pred_", p, ".RDS"))
    }

    rm(plot_list)
  }

  return(
    message(
      paste0("Finished in ", round(difftime(Sys.time(), start, unit = "secs")[[1]], digits = 1), " seconds.")
      )
    )
}
