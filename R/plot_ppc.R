#' Plot posterior predictions against observed data
#'
#' \code{plot_ppc} plots posterior predictions in a variety of ways.
#'
#' @param train_indiv List, maximum length 3. The first element should be an
#' individual-level \code{tibble} containing summed predictions for each trial
#' and individual (outputted from [get_preds_by_chain]). The second and third
#' elements should be integers or numeric vectors containing the number of
#' trials to lag for the training plots; and the last n trials to calculate
#' differences in mean observed/predicted densities for.
#' @param train_trials List, maximum length 3. The first element should be a
#' trial-level \code{tibble} containing summed posterior draws and their HDIs,
#' both overall and for each block and block group of interest (outputted
#' from [get_preds_by_chain]).
#' @param test_perf List, maximum length 3. The first element should be a
#' individual-level \code{tibble} containing summed predictions for each trial
#' and individual (outputted from [get_preds_by_chain]). The second and third
#' lists are optional, and are passed to the \code{plt.test} argument of
#' [plot_import] to plot observed grouped and individual pair accuracy
#' respectively against their posterior predictions (a grouped plot including
#' all pairs is plotted by default).
#' @param id subjID to select if only plots for a single participant are
#' desired. Will also accept a single numeric value i, which will select the
#' ith participant in the output.
#' @param group_title Sets consistent titles for all plots.
#' @param legend_pos Enables the legend positions to be set manually.
#' @param pal,font,font_size Same as [plot_import()].
#' @param ... Other rarely used arguments which set the number of trials/blocks
#' or the name of the predicted variable.
#'
#' @return Either a single or named \code{list} of \code{ggplot} objects
#'
#' @importFrom magrittr %>%
#' @importFrom rlang := !!
#' @export

plot_ppc <- function(
  train_indiv = list(),
  train_trials = list(),
  test_perf = list(),
  id = NULL,
  group_title = "",
  legend_pos = "right",
  pal = NULL,
  font = "",
  font_size = 14,
  ...) {

  if (is.null(pal)) {
    pal <- c("#ffc9b5", "#648767", "#b1ddf1", "#95a7ce", "#987284", "#3d5a80",
             "#C4F7A1", "#B1740F")
  }
  else if (!is.null(pal) & length(pal) < 8) {
    message("Need at least 8 colours, reverting to defaults.")
    pal <- c("#ffc9b5", "#648767", "#b1ddf1", "#95a7ce", "#987284", "#3d5a80",
             "#C4F7A1", "#B1740F")
  }

  ## useless assignments to appease R CMD check
  choice <- trial_no <- trial_no_group <- choice_type <- type <- acc_type <-
    subjID <- cuml_acc_mean <- cuml_acc_mean_sub_se <- cuml_acc_mean_pl_se <-
    choice_pred_prop <- obs <- post_mean_pred <- mean_obs_type <-
    mean_pred_type <- ..count.. <- avg_type <- obs_mean <- pred_post_mean <-
    pred_post_lower_95_hdi <- pred_post_upper_95_hdi <- NULL

  l <- list(...)
  if (is.null(l$max_trials_grp)) l$max_trials_grp <- 120
  if (is.null(l$block_size)) l$block_size <- 20
  if (is.null(l$out_dir)) l$out_dir <- ""
  if (is.null(l$pred_var)) l$pred_var <- "y_pred"
  if (is.null(l$n_test_trials)) l$n_test_trials <- 60

  pairs <- list("AB", "CD", "EF")
  names(pairs) <- c("12", "34", "56")
  std <- function(x) sd(x, na.rm = TRUE)/sqrt(length(!is.na(x)))

  plt_list <- list()

  if (length(train_indiv) > 0) {
    train_indiv_df <- train_indiv[[1]] %>%
      dplyr::select(-tidyselect::contains("cuml_accuracy")) %>%
      dplyr::rename(choice_obs = choice) %>%
      tidyr::pivot_longer(
        tidyselect::contains("choice"),
        names_to = "choice_type", values_to = "choice", names_prefix = "choice_"
        ) %>%
      dplyr::arrange(trial_no) %>%
      dplyr::mutate(
        acc_type = ifelse(grepl("obs", choice_type), "Observed", "Predicted")
        ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        type = paste0(pairs[[as.character(type)]], " (", acc_type, ")")
        ) %>%
      dplyr::select(-acc_type) %>%
      dplyr::group_by(subjID, type, choice_type)

    trial_lags <- tryCatch(train_indiv[[2]], error = function(e) return(list()))
    for (lag in trial_lags) {
      col_name <- rlang::sym(paste0("cuml_accuracy_l", lag))
      train_indiv_df <- train_indiv_df %>%
        dplyr::mutate(
          !!col_name := runner::runner(
            x = choice,
            f = function(x) {sum(x, na.rm = T)/sum(!is.na(x))},
            k = lag
            )
          )
    }
    cols <-
      names(train_indiv_df)[startsWith(names(train_indiv_df), "cuml_accuracy")]
    tr_plts <- list()
    for (lag_num in seq_along(trial_lags)) {
      n_lag <- trial_lags[lag_num]
      col <- rlang::sym(cols[lag_num])
      plt_name <- paste0("training_lag", n_lag)

      if (!is.null(id)) {
        train_indiv_df <- train_indiv_df %>%
          dplyr::filter(subjID == id)
      }

      tr_plot_df <- train_indiv_df %>%
        dplyr::group_by(type, trial_no_group) %>%
        dplyr::mutate(cuml_acc_mean = mean(!!col, na.rm = T)) %>%
        dplyr::mutate(cuml_acc_mean_sub_se = cuml_acc_mean - std(!!col)) %>%
        dplyr::mutate(cuml_acc_mean_pl_se = cuml_acc_mean + std(!!col)) %>%
        dplyr::ungroup() %>%
        dplyr::distinct(
          trial_no_group, type, cuml_acc_mean, cuml_acc_mean_sub_se,
          cuml_acc_mean_pl_se
        )

      plt_tr <- tr_plot_df %>%
        ggplot2::ggplot(ggplot2::aes(x = trial_no_group, y = cuml_acc_mean,
                                 colour = factor(type), fill = factor(type))) +
        ggplot2::geom_point(alpha=0.65) +
        ggplot2::geom_line() +
        ggplot2::geom_ribbon(ggplot2::aes(
          ymin = cuml_acc_mean_sub_se, ymax = cuml_acc_mean_pl_se), alpha = 0.2
        ) +
        ggplot2::scale_x_continuous(
          breaks = seq(0, l$max_trials_grp, l$block_size)
          ) +
        ggplot2::geom_vline(
          xintercept = tryCatch(c(seq(n_lag, 120 - n_lag, n_lag)),
                                error = function(e) NULL),
          linetype = "dashed", alpha = 0.5
        ) +
        ggplot2::xlab("Trial number") +
        ggplot2::ylab("Cumulative A/C/E choice probability (\u00B1 SE)") +
        ggplot2::scale_color_manual(name = "Trial Type", values = pal) +
        ggplot2::scale_fill_manual(name = "Trial Type", values = unlist(pal)) +
        cowplot::theme_half_open(
          font_size = font_size,
          font_family = font
        ) +
        ggplot2::theme(legend.position = legend_pos)

      if (n_lag == l$max_trials_grp) {
        plt_tr <- plt_tr +
          ggplot2::ggtitle(group_title, subtitle = "All trials")
      }
      else {
        plt_tr <- plt_tr +
          ggplot2::ggtitle(
            group_title, subtitle = paste0(n_lag, "-trial lagged")
            )
      }

      if (is.null(id)) {
        plt_tr <- plt_tr +
          ggplot2::geom_ribbon(ggplot2::aes(
            ymin = cuml_acc_mean_sub_se, ymax = cuml_acc_mean_pl_se),
            alpha = 0.2
          )
        }
      tr_plts[[plt_name]] <- plt_tr
    }

    if (length(tr_plts) > 0) plt_list$training_cum_prob <- tr_plts

    overall_avgs <- tryCatch(
      train_indiv[[3]], error = function(e) return(list())
      )
    if (length(overall_avgs) > 0) {
      avg_plts <- list()
      if (!is.null(id)) {
        avg_overall_df <- train_indiv[[1]] %>%
          dplyr::filter(subjID == id)
      }
      else avg_overall_df <- train_indiv[[1]]

      for (avg_diff in overall_avgs) {
        avg_overall_df <- train_indiv[[1]] %>%
          dplyr::select(subjID, trial_no_group, type, choice,
                        choice_pred_prop) %>%
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
        avg_plt <- avg_overall_df %>%
          ggplot2::ggplot(ggplot2::aes(x = diff, fill = type, colour = type)) +
          ggplot2::geom_density(ggplot2::aes(y = ..count..), alpha = 0.4) +
          ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
          ggplot2::xlab("Difference in mean A/C/E choice probability") +
          ggplot2::ylab("Count") +
          ggplot2::scale_color_manual(name = "Trial Type", values = pal) +
          ggplot2::scale_fill_manual(
            name = "Trial Type", values = unlist(pal)
            ) +
          cowplot::theme_half_open(
            font_size = font_size,
            font_family = font
          ) +
          ggplot2::theme(legend.position = legend_pos)

        if (avg_diff == l$max_trials_grp) {
          avg_plts[[avg_nm]] <- avg_plt +
            ggplot2::ggtitle(
              group_title, subtitle = "All trials (observed minus predicted)"
              )
        }
        else if (avg_diff == l$block_size) {
          avg_plts[[avg_nm]] <- avg_plt +
            ggplot2::ggtitle(
              group_title, subtitle = "Final block (observed minus predicted)"
              )
        }
        else {
          avg_plts[[avg_nm]] <- avg_plt +
            ggplot2::ggtitle(
              group_title,
              subtitle = paste0(
                "Last ", avg_diff, " trials (observed minus predicted)"
                )
              )
        }
      }
      plt_list$diffs_obs_pred <- avg_plts
    }
  }
  if (length(train_trials) > 0) {
    if (!is.null(id)) {
      avg_overall_df <- train_trials[[1]] %>%
        dplyr::filter(subjID == id)
    }
    else train_trials_df <- train_trials[[1]]

    train_trials_df <- train_trials[[1]] %>%
      dplyr::rowwise() %>%
      dplyr::mutate(avg_type = strsplit(sub("_", "\01", type),
                                        "\01")[[1]][2]) %>%
      dplyr::mutate(type = strsplit(sub("_", "\01", type), "\01")[[1]][1])

    trial_groups <- tryCatch(
      train_trials[[2]], error = function(e) return(list())
      )
    trial_plt_list <- list()
    for (trgrp in trial_groups) {
      skip_to_next <- FALSE
      tryCatch(
        plot_trials_df <- train_trials_df %>% dplyr::filter(avg_type == trgrp),
        error = function(e) skip_to_next <<- TRUE
        )
      if (!skip_to_next) {
        trial_plt_list[[trgrp]] <- plot_trials_df %>%
          ggplot2::ggplot(ggplot2::aes(x = obs_mean, y = pred_post_mean,
                                       colour = type)) +
          ggplot2::geom_point(size = 2, alpha = 0.25) +
          ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
          ggplot2::geom_errorbar(
            ggplot2::aes(x = obs_mean, ymin = pred_post_lower_95_hdi,
                         ymax = pred_post_upper_95_hdi),
            width = 0.01, alpha = 0.1
          ) +
          ggplot2::xlab("Observed mean A/C/E choice probability") +
          ggplot2::ylab(
            "Predicted mean A/C/E choice probability (\u00B1 95% HDI)"
            ) +
          ggplot2::scale_color_manual(name = "Trial Type", values = pal) +
          ggplot2::scale_fill_manual(
            name = "Trial Type", values = unlist(pal)
            ) +
          cowplot::theme_half_open(
            font_size = font_size,
            font_family = font
          ) +
          ggplot2::theme(legend.position = legend_pos) +
          ggplot2::ggtitle(
            group_title,
            subtitle = bquote(R^2~"="~.(
              round(
                cor(plot_trials_df$obs_mean, plot_trials_df$pred_post_mean)^2, 2
                )
              )~"("*.(substr(trgrp, 1, 1))*.(
                gsub("_", " ", substr(trgrp, 2, nchar(trgrp)))
                )*")"
          )
        )
      }
    }
    plt_list$indiv_posteriors <- trial_plt_list
  }
  if (length(test_perf) > 0) {

    pair_groups <- tryCatch(test_perf[[2]], error = function(e) return(list()))
    indiv_pairs <- tryCatch(test_perf[[3]], error = function(e) return(list()))

    test_perf_df <- test_perf[[1]] %>%
      dplyr::select(-tidyselect::contains("cuml_accuracy")) %>%
      dplyr::rename(choice_obs = choice) %>%
      tidyr::pivot_longer(tidyselect::contains("choice"),
                          names_to = "choice_type", values_to = "choice",
                          names_prefix = "choice_") %>%
      dplyr::arrange(trial_no) %>%
      dplyr::mutate(
        group = ifelse(grepl("obs", choice_type), "Observed", "Predicted")
        )

    if (!is.null(id)) {
      test_perf_df <- test_perf_df %>% dplyr::filter(subjID == id)
      import_single <- TRUE
    }
    else {
      import_single <- FALSE
    }

    grouped_bar_ppc <- plot_import(
        parsed_list = NULL, types = "test", plt.test = pair_groups,
        grp_compare = "group", test_df = test_perf_df,
        import_single = import_single, legend_pos = legend_pos, pal = pal,
        font = font, font_size = font_size
      ) +
      ggplot2::ggtitle(group_title, subtitle = "Test performance (grouped)")

    if (length(indiv_pairs) > 0) {
      indiv_bar_ppc <- plot_import(
        parsed_list = NULL, types = "test", plt.test = indiv_pairs,
        grp_compare = "group", test_df = test_perf_df,
        import_single = import_single, legend_pos = legend_pos, pal = pal,
        font = font, font_size = font_size
      ) +
      ggplot2::ggtitle(
        group_title, subtitle = "Test performance (individual pairs)"
        )

      plt_list$test_perf <-
        cowplot::plot_grid(
          grouped_bar_ppc + ggplot2::theme(legend.position = "none"),
          indiv_bar_ppc, nrow = 1, rel_widths = c(1,1.4)
        )
    }
    else {
      plt_list$test_perf <- grouped_bar_ppc
    }
  }
  return(plt_list)
}

