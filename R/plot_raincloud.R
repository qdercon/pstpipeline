#' Raincloud plots
#'
#' \code{plot_raincloud} plots simple "raincloud" plots of posterior densities,
#' most likely for model parameters (by group).
#'
#' @param summary_df List of [cmdstanr::summary()] outputs for the fit(s) of
#' interest.
#' @param raw_df List of raw data inputs to the above fits (in the same order).
#' Used to correctly link subject IDs to independent variables.
#' @param type Type of plot to retunn - either separate plots for each
#' \code{parameter}, or each transdiagnostic rlang::symptom \code{factor}.
#' @param test Boolean indicating whether summaries are from the test phase.
#' @param alpha_par_nms Vector of learning rate parameter names for dual
#' learning rate model. Defaults to "pos" and "neg" if not specified.
#' @param by Separately plot distributions by a certain demographic variable?
#' @param rhat_upper,ess_lower Upper and lower bounds for individuals' rhat and
#' effective sample size (ESS) for the parameters of interest. For consistency,
#' if one of an individual's parameters do not meet these criteria, all their
#' parameters are removed from the plot. Set to \code{Inf} and \code{0} to avoid
#' excluding any individuals.
#' @param legend_title,legend_labels,legend_pos Controls to name and label the
#' items in the legend (as these may be formatted poorly if \code{by != NULL}),
#' plus to set its position.
#' @param factor_scores \code{data.frame} with the derived transdiagnostic
#' factor scores, required if \code{type = "factor"}.
#' @param flip Boolean indicating if the axes should be flipped.
#' @param cred Vector, length 2, which defines the % HDI covered by the boxplot
#' boxes and lines respectively.
#' @param dist_nudge Controls the position of the distribution plots in the
#' x-direction.
#' @param box_width Control the total width of the boxplot(s).
#' @param points Either a named vector of point size, width, and nudge, or
#' \code{NULL} to points from the plot entirely.
#' @param pal,font_size,font Same as [plot_import()].
#'
#' @importFrom magrittr %>%
#' @importFrom rlang := !!
#' @importFrom stats setNames
#' @export

plot_raincloud <- function(summary_df,
                           raw_df,
                           type = "parameter",
                           test = FALSE,
                           alpha_par_nms = NULL,
                           by = NULL,
                           rhat_upper = 1.1,
                           ess_lower = 100,
                           legend_title = by,
                           legend_labels = NULL,
                           legend_pos = "right",
                           factor_scores = NULL,
                           flip = TRUE,
                           cred = c(0.95, 0.99),
                           dist_nudge = 0.1,
                           box_width = 0.6,
                           points =
                             c("size" = 0.25, "width" = 0.15, "nudge" = 0.225),
                           pal = NULL,
                           font_size = 11,
                           font = "") {

  if (is.null(pal)) {
    pal <- c("#ffc9b5", "#648767", "#b1ddf1", "#95a7ce", "#987284", "#3d5a80")
  }

  cred <- sort(cred)
  cred_l1 <- (1-cred[2])/2
  cred_l2 <- (1-cred[1])/2

  posterior_mean <- subjID <- value <- NULL # appease R CMD check

  all_data <- list()
  for (s in seq_along(summary_df)) {
    all_data[[s]] <- make_par_df(
      raw_df[[s]], summary_df[[s]], rhat_upper = rhat_upper,
      ess_lower = ess_lower
    )
  }
  all_data <- data.table::rbindlist(all_data, use.names = TRUE)
  if (!is.null(factor_scores)) {
    all_data <- all_data %>%
      dplyr::left_join(factor_scores, by = "subjID")
  }

  if (type == "parameter") {
    df <- all_data %>%
      dplyr::rename(value = posterior_mean)
    if (length(unique(df[[type]])) == 2) df[[type]] <- factor(df[[type]])
    else {
      df[[type]] <-
        factor(df[[type]], levels = c("alpha_pos", "alpha_neg", "beta"))
    }
  }
  else if (type == "factor") {
    df <- all_data %>%
      tidyr::pivot_longer(cols = c("AD", "Compul", "SW"), names_to = type) %>%
      dplyr::distinct(subjID, factor, .keep_all = T)
    if (!flip) {
      df[[type]] <- factor(df[[type]], levels = c("AD", "Compul", "SW"))
    }
    else df[[type]] <- factor(df[[type]], levels = c("SW", "Compul", "AD"))
  }

  type <- rlang::sym(type)

  if (is.null(by)) {
    rain_plot <- df %>%
      ggplot2::ggplot(ggplot2::aes(x = !!type, y = value, fill = !!type,
                                   colour = !!type)) +
      ggplot2::guides(colour = "none", fill = "none") +
      ggplot2::scale_colour_manual(values = pal) +
      ggplot2::scale_fill_manual(values = pal)
  }
  else {
    by <- rlang::sym(by)
    by_length <- length(unique(df[[by]]))
    rain_plot <- df %>%
      ggplot2::ggplot(ggplot2::aes(x = !!type, y = value, fill = factor(!!by),
                                   colour = interaction(!!by, !!type))) +
      ggplot2::guides(colour = "none") +
      ggplot2::scale_colour_manual(
        values = pal[c(1:by_length, 1:by_length, 1:by_length)]
      )
    if (!is.null(legend_labels)) {
      rain_plot <- rain_plot +
        ggplot2::scale_fill_manual(values = pal, name = legend_title,
                                   labels = legend_labels)
    } else {
      rain_plot <- rain_plot +
        ggplot2::scale_fill_manual(values = pal, name = legend_title)
    }
  }

  rain_plot <- rain_plot +
    geom_flat_violin(
      position = ggplot2::position_nudge(x = dist_nudge, y = 0),
      adjust = 2,
      trim = FALSE
    )

  if (!is.null(points)) {
    rain_plot <- rain_plot +
      ggplot2::geom_point(
        ggplot2::aes(x = as.numeric(!!type) - points[3]),
        position = ggplot2::position_jitter(width = points[2], height = 0),
        size = points[1]
      )
  }

  rain_plot <- rain_plot +
    ggplot2::stat_summary(
      geom = "boxplot",
      fun.data = function(x) {
        setNames(
          quantile_hdi(
            x, c(cred_l1, cred_l2, 0.5, 1 - cred_l2, 1 - cred_l1),
            transform = FALSE), c("ymin", "lower", "middle", "upper", "ymax")
          )
      },
      position = ggplot2::position_dodge2(),
      alpha = 0.6,
      width = box_width
    ) +
    cowplot::theme_half_open(
      font_size = font_size,
      font_family = font
    ) +
    ggplot2::theme(legend.position = legend_pos)


  if (flip) {
    rain_plot <- rain_plot +
      ggplot2::coord_flip()
  }
  if (type == "parameter") {
    if (length(unique(df[[type]])) == 2) {
      if (test) labs <- c(expression(alpha*minute), expression(beta*minute))
      else labs <- c(expression(alpha), expression(beta))
    }
    else {
      if (is.null(alpha_par_nms)) alpha_par_nms <- c("pos", "neg")
      if (test) {
        labs <-
          c(str2expression(paste0("alpha*minute[", alpha_par_nms[1], "]")),
            str2expression(paste0("alpha*minute[", alpha_par_nms[2], "]")),
            expression(beta*minute))
      }
      else {
        labs <-
          c(str2expression(paste0("alpha[", alpha_par_nms[1], "]")),
            str2expression(paste0("alpha[", alpha_par_nms[2], "]")),
            expression(beta))
      }
    }
    rain_plot <- rain_plot +
      ggplot2::scale_x_discrete(name = NULL, labels = labs) +
      ggplot2::scale_y_continuous(
        name ="Posterior mean", breaks = c(0, 0.5, 1, 2, 4, 8),
        trans = "pseudo_log"
        )
  }
  else {
    if (flip) {
      x_labels <- c("Social withdrawal", "Compulsive behaviour",
                    "Anxiety/depression")
    }
    else {
      x_labels <- c("Anxiety/depression", "Compulsive behaviour",
                    "Social withdrawal")
    }
    rain_plot <- rain_plot +
      ggplot2::scale_x_discrete(name = NULL, labels = x_labels) +
      ggplot2::ylab("Predicted factor score")
  }
  return(rain_plot)
}
