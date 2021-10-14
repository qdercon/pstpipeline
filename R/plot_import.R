#' Plot raw experiment data
#'
#' \code{plot_import} plots a single participant's data from \code{pstpipeline::import_single()}, or
#' a single participant (if \code{!is.null(id)}) or all participants' data from
#' \code{pstpipeline::import_multiple()}.
#'
#' @param parsed_df \code{pstpipeline::import_single()} or \code{pstpipeline::import_multiple()} output.
#' @param import_single Is the output from \code{pstpipeline::import_single()}?
#' @param id Prolific or other participant ID to select if only a single partipant's data is desired.
#' @param pal Define a custom colour palette for the plots? Otherwise reverts to defaults.
#' @param font Use a custom font for the plots? Will likely require \code{extrafont::font_import()} to
#' be run first.
#' @param font_size Base plot font size passed to \code{ggplot2::theme_gray}.
#' @param plot_type Select plots as required.
#'
#' @return \code{list} of \code{ggplot} objects
#'
#' @importFrom magrittr %>%
#' @export

plot_import <-
  function(parsed_df, import_single = FALSE, id = NULL, pal = NULL, font = "", font_size = 14,
           plot_type = c("tr20", "tr60", "tr_all", "tr_questions", "happy", "confident","engaged",
                         "test_perf")) {

    if (is.null(pal)) {
      pal <- c("#ffc9b5", "#648767", "#b1ddf1", "#95a7ce", "#987284", "#3d5a80")
    } else if (!is.null(pal) & length(pal) < 4) {
      message("Need at least 4 colours, reverting to defaults.")
      pal <- c("#ffc9b5", "#648767", "#b1ddf1", "#95a7ce", "#987284", "#3d5a80")
    }

    if (!import_single & !is.null(id)) {
      if (is.null(parsed_df$individual_results)) {
        stop(paste0("Need individual results list to plot a single participant's result - ",
                    "please re-run parse_multiple with indiv = TRUE."))
      }

      if (is.numeric(subjID)) id <- names(parsed_df$individual_results)[[subjID]]
      else id <- paste0("ID", as.character(subjID))
      message(paste0("Plotting data for subject ", gsub("ID", "", id), "..."))

      training <- parsed_df$individual_results[[id]]$training
      test <- parsed_df$individual_results[[id]]$test
    }
    else {
      training <- parsed_df$training
      test <- parsed_df$test
    }

    ret <- list()

    if (!is.null(font)) {
      extrafont::loadfonts(device = "win", quiet = TRUE)
    }

    if (any(plot_type=="tr20")) {
      plot20 <- training %>%
        ggplot2::ggplot(ggplot2::aes(x=trial_no, y=cum_prob_l20, color=factor(type))) +
        ggplot2::geom_point(alpha=0.65) +
        ggplot2::geom_line() +
        ggplot2::geom_vline(xintercept=60, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=120, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=180, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=240, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=300, linetype="dashed", alpha=0.5) +
        ggplot2::scale_x_continuous(breaks=seq(0,360,30)) +
        ggplot2::xlab("Trial index") +
        ggplot2::ylab("Cumulative probability of picking stimulus A/C/E") +
        ggplot2::scale_color_manual(name = "Trial Type", labels = c("AB", "CD", "EF"), values = pal) +
        cowplot::theme_half_open(
          font_size = font_size,
          font_family = font
        ) +
        ggplot2::theme(legend.position = c(0.85, 0.2)) +
        ggplot2::ggtitle("20-trial lagged cumulative probabilities of picking correct stimulus")

      ret$training_lag_20 <- plot20

    }

    if (any(plot_type=="tr60")) {
      plot60 <- training %>%
        ggplot2::ggplot(ggplot2::aes(x=trial_no, y=cum_prob_l60, color=factor(type))) +
        ggplot2::geom_point(alpha=0.65) +
        ggplot2::geom_line() +
        ggplot2::geom_vline(xintercept=60, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=120, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=180, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=240, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=300, linetype="dashed", alpha=0.5) +
        ggplot2::scale_x_continuous(breaks=seq(0,360,30)) +
        ggplot2::xlab("Trial index") +
        ggplot2::ylab("Cumulative probability of picking stimulus A/C/E") +
        ggplot2::scale_color_manual(name = "Trial Type", labels = c("AB", "CD", "EF"), values = pal) +
        cowplot::theme_half_open(
          font_size = font_size,
          font_family = font
        ) +
        ggplot2::theme(legend.position = c(0.85, 0.2)) +
        ggplot2::ggtitle("60-trial lagged cumulative probabilities of picking correct stimulus")

      ret$training_lag_60 <- plot60

    }

    if (any(plot_type=="tr_all")) {
      plot_full <- training %>%
        ggplot2::ggplot(ggplot2::aes(x=trial_no, y=cum_prob_all, color=factor(type))) +
        ggplot2::geom_point(alpha=0.65) +
        ggplot2::geom_line() +
        ggplot2::geom_vline(xintercept=60, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=120, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=180, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=240, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=300, linetype="dashed", alpha=0.5) +
        ggplot2::scale_x_continuous(breaks=seq(0,360,30)) +
        ggplot2::xlab("Trial index") +
        ggplot2::ylab("Cumulative probability of picking stimulus A/C/E") +
        ggplot2::scale_color_manual(name = "Trial Type", labels = c("AB", "CD", "EF"), values = pal) +
        cowplot::theme_half_open(
          font_size = font_size,
          font_family = font
        ) +
        ggplot2::theme(legend.position = c(0.85, 0.2)) +
        ggplot2::ggtitle("Overall cumulative probabilities of picking correct stimulus")

      ret$training_all <- plot_full

    }

    if (any(plot_type=="tr_questions")) {
      fatigue_questions <- training %>%
        dplyr::select(subjID, trial_block, trial_no, fatigue_rt, fatigue_slider_start,
                      fatigue_response) %>%
        dplyr::rename(question_rt = fatigue_rt, question_slider_start = fatigue_slider_start,
              question_response = fatigue_response) %>%
        dplyr::mutate(question_type="fatigue") %>%
        tidyr::drop_na()

      training_questions <- training %>%
        dplyr::select(subjID, trial_block, trial_no, question_type, question_slider_start, question_rt,
                  question_response)

      plot_tr_q <- fatigue_questions %>%
        dplyr::bind_rows(training_questions) %>%
        ggplot2::ggplot(
          ggplot2::aes(
            x=trial_no, y=question_response,
            color = factor(question_type, levels=c("happy", "confident", "engaged", "fatigue")))
          ) +
        ggplot2::geom_point(alpha=0.65) +
        ggplot2::geom_line() +
        ggplot2::geom_vline(xintercept=60, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=120, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=180, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=240, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=300, linetype="dashed", alpha=0.5) +
        ggplot2::scale_x_continuous(breaks=seq(0,360,30)) +
        ggplot2::xlab("Trial index") +
        ggplot2::ylab("Rating (/100)") +
        ggplot2::scale_color_manual(
          name = "Affect noun",
          labels = c("Happiness", "Confidence", "Engagement", "Fatigue"),
          values = pal
        ) +
        cowplot::theme_half_open(
          font_size = font_size,
          font_family = font
        ) +
        ggplot2::ggtitle("Subjective affect over training")

      ret$affect_questions <- plot_tr_q
    }

    if (any(plot_type=="happy")) {
      plot_happy <- training %>%
        dplyr::filter(question_type=="happy") %>%
        ggplot2::ggplot(ggplot2::aes(x=trial_no, y=question_response,
                                     color=factor(correct, levels=c("TRUE", "FALSE")))) +
        ggplot2::geom_point(alpha=0.65) +
        ggplot2::geom_line() +
        ggplot2::geom_vline(xintercept=60, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=120, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=180, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=240, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=300, linetype="dashed", alpha=0.5) +
        ggplot2::scale_x_continuous(breaks=seq(0,360,30)) +
        ggplot2::xlab("Trial index") +
        ggplot2::ylab("Rating (/100)") +
        ggplot2::scale_color_manual(name = "Rewarded?", labels = c("Yes", "No"), values = pal) +
        cowplot::theme_half_open(
          font_size = font_size,
          font_family = font
        ) +
        ggplot2::ggtitle("Subjective happiness in rewarded and non-rewarded trials")

      ret$happy <- plot_happy
    }

    if (any(plot_type=="engaged")) {
      plot_engaged <- training %>%
        dplyr::filter(question_type=="engaged") %>%
        ggplot2::ggplot(ggplot2::aes(x=trial_no, y=question_response,
                                     color=factor(correct, levels=c("TRUE", "FALSE")))) +
        ggplot2::geom_point(alpha=0.65) +
        ggplot2::geom_line() +
        ggplot2::geom_vline(xintercept=60, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=120, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=180, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=240, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=300, linetype="dashed", alpha=0.5) +
        ggplot2::scale_x_continuous(breaks=seq(0,360,30)) +
        ggplot2::xlab("Trial index") +
        ggplot2::ylab("Rating (/100)") +
        ggplot2::scale_color_manual(name = "Rewarded?", labels = c("Yes", "No"), values = pal) +
        cowplot::theme_half_open(
          font_size = font_size,
          font_family = font
        ) +
        ggplot2::ggtitle("Subjective engagement in rewarded and non-rewarded trials")

      ret$engaged <- plot_engaged
    }

    if (any(plot_type=="confident")) {
      plot_conf <- training %>%
        dplyr::filter(question_type=="confident") %>%
        ggplot2::ggplot(ggplot2::aes(x=trial_no, y=question_response,
                                     color=factor(correct, levels=c("TRUE", "FALSE")))) +
        ggplot2::geom_point(alpha=0.65) +
        ggplot2::geom_line() +
        ggplot2::geom_vline(xintercept=60, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=120, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=180, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=240, linetype="dashed", alpha=0.5) +
        ggplot2::geom_vline(xintercept=300, linetype="dashed", alpha=0.5) +
        ggplot2::scale_x_continuous(breaks=seq(0,360,30)) +
        ggplot2::xlab("Trial index") +
        ggplot2::ylab("Rating (/100)") +
        ggplot2::scale_color_manual(name = "Rewarded?", labels = c("Yes", "No"), values = pal) +
        cowplot::theme_half_open(
          font_size = font_size,
          font_family = font
        ) +
        ggplot2::ggtitle("Subjective confidence in rewarded and non-rewarded trials")

      ret$confidence <- plot_conf
    }

    if (any(plot_type=="test_perf")) {
      plot_test <- test %>%
        ggplot2::ggplot(ggplot2::aes(x=factor(test_type,
                                              levels=c("chooseA","avoidB", "novel", "training")),
                   fill=factor(correct, levels=c("TRUE", "FALSE")))) +
        ggplot2::geom_bar() +
        ggplot2::geom_text(stat = "count", family = ifelse(!is.null(font), font, ""),
                           ggplot2::aes(label = ggplot2::after_stat(count),
                                        colour=factor(correct, levels=c("TRUE", "FALSE"))),
                  position = ggplot2::position_stack(vjust=0.5)) +
        ggplot2::xlab("Test type") +
        ggplot2::ylab("Count") +
        ggplot2::scale_x_discrete() +
        ggplot2::scale_fill_manual(values = pal, name = NULL, labels=c("Correct", "Incorrect")) +
        ggplot2::scale_colour_manual(values=c("#000000", "#FFFFFF"), guide= "none") +
        ggplot2::scale_alpha_manual(values=0.85) +
        cowplot::theme_half_open(
          font_size = font_size,
          font_family = font
        ) +
        ggplot2::ggtitle("Test performance")

      ret$testperf <- plot_test
    }

  return(ret)
}
