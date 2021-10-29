#' Plot posterior predictions against observed data
#'
#' \code{plot_ppc} plots
#'
#' @param train_indiv List, maximum length 3. The first element should be an individual-level \code{tibble}
#' containing summed predictions for each trial and individual (outputted from [pstpipeline::get_preds_by_chain]).
#' The second and third elements should be integers or numeric vectors containing the number of trials to lag for
#' the training plots; and the last n trials to calculate differences in mean observed/predicted densities for.
#' @param train_trials List, maximum length 3. The first element should be a trial-level \code{tibble} containing
#' summed posterior draws and their HDIs, both overall and for each block and block group of interest (outputted
#' from [pstpipeline::get_preds_by_chain]).
#' @param id subjID to select if only plots for a single participant are desired. Will also accept a single
#' numeric value i, which will select the ith participant in the output.
#' @param legend_pos Enables the legend positions to be set manually.
#' @param pal Define a custom colour palette for the plots? Otherwise reverts to defaults.
#' @param font Use a custom font for the plots? Warnings suggest \code{extrafont::font_import()} should be run.
#' @param font_size Base plot font size.
#'
#' @return Either a single or named \code{list} of \code{ggplot} objects
#'
#' @importFrom magrittr %>%
#' @importFrom rlang := !!
#' @export

plot_ppc <- function(
  train_indiv = list(),
  train_trials = list(),
  test_phase = list(),
  id = NULL,
  legend_pos = "right",
  pal = NULL,
  font = "",
  font_size = 14,
  ...) {

  if (is.null(pal)) {
    pal <- c("#ffc9b5", "#648767", "#b1ddf1", "#95a7ce", "#987284", "#3d5a80")
  }
  else if (!is.null(pal) & length(pal) < 6) {
    message("Need at least 6 colours, reverting to defaults.")
    pal <- c("#ffc9b5", "#648767", "#b1ddf1", "#95a7ce", "#987284", "#3d5a80")
  }
  if (!is.null(font)) {
    extrafont::loadfonts(device = "win", quiet = TRUE)
  }

  l <- list(...)
  if (is.null(l$max_trials_grp)) l$max_trials_grp <- 120
  if (is.null(l$block_size)) l$block_size <- 20
  if (is.null(l$out_dir)) l$out_dir <- ""
  if (is.null(l$pred_var)) l$pred_var <- "y_pred"

  pairs <- list("AB", "CD", "EF")
  names(pairs) <- c("12", "34", "56")
  std <- function(x) sd(x, na.rm = TRUE)/sqrt(length(!is.na(x)))

  plt_list <- list()

  if (length(train_indiv) > 0) {
    train_indiv_df <- train_indiv[[1]] %>%
      dplyr::select(-contains("cuml_accuracy")) %>%
      dplyr::rename(choice_obs = choice) %>%
      tidyr::pivot_longer(contains("choice"), names_to = "choice_type", values_to = "choice",
                          names_prefix = "choice_") %>%
      dplyr::arrange(trial_no) %>%
      dplyr::mutate(acc_type = ifelse(grepl("obs", choice_type), "Observed", "Predicted")) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(type = paste0(pairs[[as.character(type)]], " (", acc_type, ")")) %>%
      dplyr::select(-acc_type) %>%
      dplyr::group_by(subjID, type, choice_type)

    trial_lags <- tryCatch(train_indiv[[2]], error = function(e) return(list()))
    for (lag in trial_lags) {
      col_name <- rlang::sym(paste0("cuml_accuracy_l", lag))
      train_indiv_df <- train_indiv_df %>%
        dplyr::mutate(
          !!col_name := runner::runner(
            x = choice, f = function(x) {sum(x, na.rm = T)/sum(!is.na(x))}, k = lag)
          )
    }
    cols <- names(train_indiv_df)[startsWith(names(train_indiv_df), "cuml_accuracy")]
    tr_plts <- list()
    for (lag_num in seq_along(trial_lags)) {
      n_lag <- trial_lags[lag_num]
      col <- rlang::sym(cols[lag_num])
      plt_name <- paste0("training_lag", n_lag)

      if (!is.null(id)) train_indiv_df <- train_indiv_df %>% dplyr::filter(subjID == id)

      tr_plot_df <- train_indiv_df %>%
        dplyr::group_by(type, trial_no_group) %>%
        dplyr::mutate(cuml_acc_mean = mean(!!col, na.rm = T)) %>%
        dplyr::mutate(cuml_acc_mean_sub_se = cuml_acc_mean - std(!!col)) %>%
        dplyr::mutate(cuml_acc_mean_pl_se = cuml_acc_mean + std(!!col)) %>%
        dplyr::ungroup() %>%
        dplyr::distinct(
          trial_no_group, type, cuml_acc_mean, cuml_acc_mean_sub_se, cuml_acc_mean_pl_se
        )

      plt_tr <- tr_plot_df %>%
        ggplot2::ggplot(ggplot2::aes(x = trial_no_group, y = cuml_acc_mean,
                                 colour = factor(type), fill = factor(type))) +
        ggplot2::geom_point(alpha=0.65) +
        ggplot2::geom_line() +
        ggplot2::geom_ribbon(ggplot2::aes(
          ymin = cuml_acc_mean_sub_se, ymax = cuml_acc_mean_pl_se), alpha = 0.2
        ) +
        ggplot2::scale_x_continuous(breaks = seq(0, l$max_trials_grp, l$block_size)) +
        ggplot2::geom_vline(
          xintercept = tryCatch(c(seq(n_lag, 120 - n_lag, n_lag)), error = function(e) NULL),
          linetype = "dashed", alpha = 0.5
        ) +
        ggplot2::xlab("Trial number") +
        ggplot2::ylab("Cumulative A/C/E choice probability (± SE)") +
        ggplot2::scale_color_manual(name = "Trial Type", values = pal) +
        ggplot2::scale_fill_manual(name = "Trial Type", values = unlist(pal)) +
        cowplot::theme_half_open(
          font_size = font_size,
          font_family = font
        ) +
        ggplot2::theme(legend.position = legend_pos) +
        ggplot2::ggtitle(paste0(n_lag, "-trial lagged"))

      if (is.null(id)) {
        plt_tr <- plt_tr +
          ggplot2::geom_ribbon(ggplot2::aes(
            ymin = cuml_acc_mean_sub_se, ymax = cuml_acc_mean_pl_se), alpha = 0.2
          )
        }
      tr_plts[[plt_name]] <- plt_tr
    }

    if (length(tr_plts) > 0) plt_list$training_cum_prob <- tr_plts

    overall_avgs <- tryCatch(train_indiv[[3]], error = function(e) return(list()))
    if (length(overall_avgs) > 0) {
      avg_plts <- list()
      if (!is.null(id)) avg_overall_df <- train_indiv[[1]] %>% dplyr::filter(subjID == id)
      else avg_overall_df <- train_indiv[[1]]

      for (avg_diff in overall_avgs) {
        avg_overall_df <- train_indiv[[1]] %>%
          dplyr::select(subjID, trial_no_group, type, choice, choice_pred_prop) %>%
          dplyr::rename(obs = choice, post_mean_pred = choice_pred_prop) %>%
          dplyr::filter(trial_no_group >= (l$max_trials_grp - avg_diff)) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(type = pairs[[as.character(type)]]) %>%
          dplyr::group_by(subjID, type) %>%
          dplyr::mutate(mean_obs_type = mean(obs)) %>%
          dplyr::mutate(mean_pred_type = mean(post_mean_pred)) %>%
          dplyr::mutate(diff = mean_obs_type - mean_pred_type) %>%
          dplyr::distinct(type, diff)

        avg_nm <- paste("last", avg_diff, "trials", sep = "_")
        avg_plts[[avg_nm]] <- avg_overall_df %>%
          ggplot2::ggplot(ggplot2::aes(x = diff, fill = type, colour = type)) +
          ggplot2::geom_density(ggplot2::aes(y = ..count..), alpha = 0.4) +
          ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
          ggplot2::xlab("Difference in mean A/C/E choice probability") +
          ggplot2::ylab("Count") +
          ggplot2::scale_color_manual(name = "Trial Type", values = pal) +
          ggplot2::scale_fill_manual(name = "Trial Type", values = unlist(pal)) +
          cowplot::theme_half_open(
            font_size = font_size,
            font_family = font
          ) +
          ggplot2::theme(legend.position = legend_pos) +
          ggplot2::ggtitle(
            paste0("Last ", avg_diff, " trials (observed minus predicted)")
            )

      }
      plt_list$diffs_obs_pred <- avg_plts
    }
  }
  if (length(train_trials) > 0) {
    if (!is.null(id)) avg_overall_df <- train_trials[[1]] %>% dplyr::filter(subjID == id)
    else train_trials_df <- train_trials[[1]]

    train_trials_df <- train_trials[[1]] %>%
      dplyr::rowwise() %>%
      dplyr::mutate(avg_type = strsplit(sub("_", "\01", type),
                                        "\01")[[1]][2]) %>%
      dplyr::mutate(type = strsplit(sub("_", "\01", type), "\01")[[1]][1])

    trial_groups <- tryCatch(train_trials[[2]], error = function(e) return(list()))
    trial_plt_list <- list()
    for (trgrp in trial_groups) {
      skip_to_next <- FALSE
      tryCatch(plot_trials_df <- train_trials_df %>% dplyr::filter(avg_type == trgrp),
               error = function(e) skip_to_next <<- TRUE)
      if (!skip_to_next) {
        trial_plt_list[[trgrp]] <- plot_trials_df %>%
          ggplot2::ggplot(ggplot2::aes(x = obs_mean, y = pred_post_mean,
                                       colour = type)) +
          ggplot2::geom_point(size = 2, alpha = 0.75) +
          ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
          ggplot2::geom_errorbar(
            ggplot2::aes(x = obs_mean, ymin = pred_post_lower_95_hdi, ymax = pred_post_upper_95_hdi),
            width = 0.01, alpha = 0.5
          ) +
          ggplot2::xlab("Observed mean A/C/E choice probability") +
          ggplot2::ylab("Predicted mean A/C/E choice probability (± 95% HDI)") +
          ggplot2::scale_color_manual(name = "Trial Type", values = pal) +
          ggplot2::scale_fill_manual(name = "Trial Type", values = unlist(pal)) +
          cowplot::theme_half_open(
            font_size = font_size,
            font_family = font
          ) +
          ggplot2::theme(legend.position = legend_pos) +
          ggplot2::ggtitle(
            paste0(toupper(substr(trgrp, 1, 1)), gsub("_", " ", substr(trgrp, 2, nchar(trgrp)))),
            subtitle = bquote(R^2~"="~.(
              round(cor(plot_trials_df$obs_mean, plot_trials_df$pred_post_mean)^2, 2)
              )
            )
          )
      }
    }

    plt_list$indiv_posteriors <- trial_plt_list

  }
  if (length(test_phase) > 0) {message("to add..")}

  return(plt_list)
}
