#' Compare learning parameters between groups over different numbers of training blocks
#'
#' @param all_res List output from [import_multiple()], split by variable of interest.
#' @param var_of_interest,covariates Same as [parameter_glm()].
#' @param model,vb,out_dir Same as [fit_learning_model()].
#' @param fit_together Whether or not to fit all participant data to the model at the
#' same time, or by group (if \code{FALSE}). Recommended for variational fits.
#' @param iter_warmup_glm,iter_sampling_glm Number of warm-up and sampling iterations
#' to fit Bayesian GLMs with (passed to [parameter_glm()]).
#' @param min_blocks,max_blocks Minimum and maximum number of blocks to fit models on.
#' For example, if \code{min_blocks} is set to 3 and \code{max_blocks} is set to 4, then
#' only fits to blocks 1 to 3 and 1 to 4 will be outputted.
#' @param ... Other arguments to pass to [fit_learning_model()] and/or [parameter_glm()].
#'
#' @importFrom magrittr %>%
#' @export

compare_block_diffs <- function(all_res,
                                var_of_interest,
                                covariates,
                                model,
                                vb = TRUE,
                                fit_together = vb,
                                out_dir = "outputs/cmdstan/compare_blocks",
                                iter_warmup_glm = 2000,
                                iter_sampling_glm = 4000,
                                min_blocks = 1,
                                max_blocks = 6,
                                ...) {

  l <- list(...)
  if (is.null(l$task_excl)) l$task_excl <- TRUE
  if (is.null(l$accuracy_excl)) l$accuracy_excl <- FALSE

  block_group <- trial_no <- NULL # to appease R CMD check

  rel_data <- list()
  if (fit_together) {
    rel_data$ppt_info <- dplyr::bind_rows(
      all_res[[1]]$ppt_info, all_res[[2]]$ppt_info
    )
    rel_data$training <- dplyr::bind_rows(
      all_res[[1]]$training, all_res[[2]]$training
    )
  } else {
    rel_data <- all_res
  }

  iter <- min_blocks:max_blocks

  par_df_ls <- list()

  for (i in seq_along(iter)) {
    if (fit_together) {
      rel_data_tr <- rel_data
      rel_data_tr$training <- rel_data_tr$training %>%
        dplyr::filter(trial_no <= iter[i]*60)

      fit_typ <- ifelse(vb, "vb", "mcmc")

      first_n_blks <- fit_learning_model(
        rel_data_tr, model = model, exp_part = "training", vb = vb,
        out_dir = out_dir, ppc = FALSE, task_excl = l$task_excl,
        accuracy_excl = l$accuracy_excl, model_checks = FALSE,
        save_model_as = paste(
          "first", i, "training_blocks", model, fit_typ, sep = "_"
          ),
        outputs = c("raw_df", "summary"), save_outputs = FALSE, ...
      )

      par_df_ls[[i]] <- parameter_glm(
        summary_df = list(first_n_blks$summary),
        raw_df = list(first_n_blks$raw_df),
        var_of_interest = var_of_interest,
        covariates = covariates,
        iter_warmup = iter_warmup_glm,
        iter_sampling = iter_sampling_glm,
        ...
      )
    }
    else {
      rel_data_gr1 <- rel_data[[1]]
      rel_data_gr2 <- rel_data[[2]]

      rel_data_gr1$training <- rel_data_gr1$training %>%
        dplyr::filter(trial_no <= iter[i]*60)
      rel_data_gr1$training <- rel_data_gr1$training %>%
        dplyr::filter(trial_no <= iter[i]*60)
      fit_typ <- ifelse(vb, "vb", "mcmc")
      grp_names <- names(rel_data)

      first_n_blks_gr1 <- fit_learning_model(
        rel_data_gr1, model = model, exp_part = "training", vb = vb,
        out_dir = out_dir, ppc = FALSE, task_excl = l$task_excl,
        accuracy_excl = l$accuracy_excl, model_checks = FALSE,
        save_model_as = paste(
          "first", i, "training_blocks", model, grp_names[1], fit_typ, sep = "_"
          ),
        outputs = c("raw_df", "summary"), save_outputs = FALSE, ...
      )
      first_n_blks_gr2 <- fit_learning_model(
        rel_data_gr2, model = model, exp_part = "training", vb = vb,
        out_dir = out_dir, ppc = FALSE, task_excl = l$task_excl,
        accuracy_excl = l$accuracy_excl, model_checks = FALSE,
        save_model_as = paste(
          "first", i, "training_blocks", model, grp_names[2], fit_typ, sep = "_"
          ),
        outputs = c("raw_df", "summary"), save_outputs = FALSE, ...
      )

      par_df_ls[[i]] <- parameter_glm(
        summary_df = list(first_n_blks_gr1$summary, first_n_blks_gr2$summary),
        raw_df = list(first_n_blk_gr1$raw_df, first_n_blk_gr2$raw_df),
        var_of_interest = var_of_interest,
        covariates = covariates,
        iter_warmup = iter_warmup_glm,
        iter_sampling = iter_sampling_glm,
        ...
      )
    }
  }

  names_all <- c("Block 1 only", sapply(2:5, function(b) paste0("Block 1 to ", b)), "All blocks")
  names(par_df_ls) <- names_all[c(min_blocks:max_blocks)]

  glm_pars_df <- data.table::rbindlist(par_df_ls, idcol = "block_group") %>%
    dplyr::mutate(block_group = factor(block_group, levels = rev(names(par_df_ls))))

  return(glm_pars_df)
}
