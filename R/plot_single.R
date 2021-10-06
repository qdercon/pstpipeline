#' Plot raw experiment data (single result)
#'
#' \code{plot_single} plots a single participant's data from \code{pstpipeline::import_multiple()}
#' output.
#'
#' @param subjID Prolific or other participant ID to select.
#' @param parsed_multiple \code{pstpipeline::import_multiple()} output
#' @param plot_type Select plots as required.
#'
#' @return \code{list} of \code{ggplot} objects
#'
#' @importFrom magrittr %>%
#' @export

plot_single <-
  function(subjID, parsed_multiple, plot_type=c("tr20", "tr60", "tr_all", "tr_questions",
                                                "happy", "confident","engaged", "test_perf")) {
    if (is.null(parsed_multiple$individual_results)) {
      stop("Need individual results list for plotting - please re-run parse_multiple with indiv=T.")
    }

    if (is.numeric(subjID)) id <- names(parsed_multiple$individual_results)[[subjID]]
    else id <- paste0("ID", as.character(subjID))
    message(paste0("Plotting data for subject ", gsub("ID", "", id), "..."))

    training <- parsed_multiple$individual_results[[id]]$training
    test <- parsed_multiple$individual_results[[id]]$test

    ret <- list()

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
        ggsci::scale_color_jama(name = "Trial Type", labels = c("AB", "CD", "EF")) +
        ggplot2::theme_grey(base_size = 14) +
        ggplot2::theme(legend.position = c(0.9, 0.1)) +
        ggplot2::ggtitle("20-trial lagged cumulative probabilities of picking correct stimulus")

      ret$plot_20 <- plot20

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
        ggsci::scale_color_jama(name = "Trial Type", labels = c("AB", "CD", "EF")) +
        ggplot2::theme_grey(base_size = 14) +
        ggplot2::theme(legend.position = c(0.9, 0.1)) +
        ggplot2::ggtitle("60-trial lagged cumulative probabilities of picking correct stimulus")

      ret$plot_60 <- plot60

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
        ggsci::scale_color_jama(name = "Trial Type", labels = c("AB", "CD", "EF")) +
        ggplot2::theme_grey(base_size = 14) +
        ggplot2::theme(legend.position = c(0.9, 0.1)) +
        ggplot2::ggtitle("Overall cumulative probabilities of picking correct stimulus")

      ret$plot_all <- plot_full

    }

    if (any(plot_type=="tr_questions")) {

      fatigue_questions <- training %>%
        dplyr::select(subjID, trial_block, trial_no, fatigue_rt, fatigue_slider_start,
                      fatigue_response) %>%
        tidyr::drop_na()
      training_questions <- training %>%
        dplyr::select(subjID, trial_block, trial_no, question_type, question_slider_start,
                      question_rt, question_response)

      plot_tr_q <- fatigue_questions %>%
        dplyr::rename(question_rt = fatigue_rt, question_slider_start = fatigue_slider_start,
                      question_response = fatigue_response) %>%
        dplyr::mutate(question_type="fatigue") %>%
        dplyr::bind_rows(training_questions) %>%
        ggplot2::ggplot(ggplot2::aes(x=trial_no, y=question_response,
                            color=factor(question_type,
                                         levels=c("happy", "confident", "engaged", "fatigue")))) +
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
        ggsci::scale_color_jama(name = "Affect noun", labels = c("Happiness", "Confidence",
                                                                 "Engagement", "Fatigue")) +
        ggplot2::theme_grey(base_size = 14) +
        ggplot2::ggtitle("Subjective affect over training")

      ret$plot_questions <- plot_tr_q
    }

    if (any(plot_type=="happy")) {
      plot_happy <- training %>%
        dplyr::filter(question_type=="happy") %>%
        ggplot2::ggplot(ggplot2::aes(x=trial_no, y=question_response, color=factor(correct,
                                                                 levels = c("TRUE", "FALSE")))) +
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
        ggsci::scale_color_jama(name = "Rewarded?", labels = c("Yes", "No")) +
        ggplot2::theme_grey(base_size = 14) +
        ggplot2::ggtitle("Subjective happiness in rewarded and non-rewarded trials")

      ret$plot_happy <- plot_happy
    }

    if (any(plot_type=="engaged")) {
      plot_engaged <- training %>%
        dplyr::filter(question_type=="engaged") %>%
        ggplot2::ggplot(ggplot2::aes(x=trial_no, y=question_response,
                            color = factor(correct, levels=c("TRUE", "FALSE")))) +
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
        ggsci::scale_color_jama(name = "Rewarded?", labels = c("Yes", "No")) +
        ggplot2::theme_grey(base_size = 14) +
        ggplot2::ggtitle("Subjective engagement in rewarded and non-rewarded trials")

      ret$plot_engaged <- plot_engaged
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
        ggsci::scale_color_jama(name = "Rewarded?", labels = c("Yes", "No")) +
        ggplot2::theme_grey(base_size = 14) +
        ggplot2::ggtitle("Subjective confidence in rewarded and non-rewarded trials")

      ret$plot_confident <- plot_conf
    }

    if (any(plot_type=="test_perf")) {
      plot_test <- test %>%
        ggplot2::ggplot(ggplot2::aes(x=factor(test_type, levels=c("chooseA","avoidB", "novel", "training")),
                            fill=factor(correct, levels=c("TRUE", "FALSE")))) +
        ggplot2::geom_bar() +
        ggplot2::geom_text(stat = "count", ggplot2::aes(label = after_stat(count),
                                      colour=factor(correct, levels=c("TRUE", "FALSE"))),
                  position = ggplot2::position_stack(vjust=0.5)) +
        ggplot2::xlab("Test type") +
        ggplot2::ylab("Count") +
        ggplot2::scale_x_discrete() +
        ggplot2::scale_fill_manual(values=ggsci::pal_jama("default")(4)[3:4], name=NULL,
                                   labels=c("Correct", "Incorrect")) +
        ggplot2::scale_colour_manual(values=c("#000000", "#FFFFFF"), guide=F) +
        ggplot2::scale_alpha_manual(values=0.85) +
        ggplot2::theme_grey(base_size = 14) +
        ggplot2::ggtitle("Test performance")

      ret$plot_testperf <- plot_test
    }

  return(ret)
}
