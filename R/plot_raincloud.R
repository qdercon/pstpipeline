#' Raincloud plots
#'
#' \code{plot_raincloud} plots simple "raincloud" plots of posterior densities,
#' most likely for model parameters (by group).
#'
#' @param summary_df List of [cmdstanr::summary()] outputs for the fit(s) of
#' interest.
#' @param raw_df List of raw data inputs to the above fits (in the same order).
#' Used to correctly link subject IDs to independent variables.
#' @param pars_of_interest Vector of parameters to plot, for example if we want
#' to only plot a specific parameter by a certain demographic variable.
#' @param type Type of plot to return - either separate plots for each
#' \code{parameter}, or each transdiagnostic symptom \code{factor}.
#' @param test Boolean indicating whether summaries are from the test phase.
#' @param by Separately plot distributions by a certain demographic variable?
#' @param recode_na Value to recode NAs in the \code{by} variable to.
#' @param par_nms The name of the column that contains the discrete axis labels
#' (i.e., the parameter or factor score names).
#' @param par_vals The name of the column that contains the continuous values to
#' plot (i.e., the parameter values or factor scores).
#' @param par_lvls The order in which to plot the parameters or factors.
#' @param par_labs Modify the labels of the parameters or factors - must be in
#' the same order as \code{par_lvls}.
#' @param afct_filt Affect adjective to filter on, if \code{affect = TRUE}.
#' @param afct_col Column containing the affect types.
#' @param combine_with A column to combine with the \code{par_nms} variable, if
#' this is further split by another variable.
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
#' @param scatter Either a named vector of scatterplot width, and nudge. Use
#' \code{NULL} to omit the scatterplots from the plot entirely.
#' @param pal,font_size,font Same as [plot_import()].
#' @param ... Additional arguments such as \code{alpha_par_nms} to alter axis
#' titles or control the \code{moment} used in [quantile_hdi()].
#'
#' @returns A ggplot.
#'
#' @examples \dontrun{
#' data(example_data)
#'
#' fit_nd <- fit_learning_model(
#'   example_data$nd,
#'   model = "2a",
#'   vb = FALSE,
#'   exp_part = "training"
#' )
#' fit_dis <- fit_learning_model(
#'   example_data$dis,
#'   model = "2a",
#'   vb = FALSE,
#'   exp_part = "training"
#' )
#'
#' plot_raincloud(
#'   summary_df = list(fit_nd$summary, fit_dis$summary),
#'   raw_df = list(fit_nd$raw_df, fit_dis$raw_df),
#'   by = "distanced",
#'   flip = FALSE
#' )
#' }
#'
#' @importFrom rlang := !!
#' @importFrom stats setNames
#' @export

plot_raincloud <- function(summary_df,
                           raw_df,
                           pars_of_interest = NULL,
                           type = "parameter",
                           test = FALSE,
                           by = NULL,
                           recode_na = 0,
                           par_nms = type,
                           par_vals = "posterior_mean",
                           par_lvls = NULL,
                           par_labs = NULL,
                           afct_filt = NULL,
                           afct_col = NULL,
                           combine_with = NULL,
                           legend_title = by,
                           legend_labels = NULL,
                           legend_pos = "right",
                           factor_scores = NULL,
                           flip = TRUE,
                           cred = c(0.95, 0.99),
                           dist_nudge = 0.1,
                           box_width = 0.1,
                           scatter = c("width" = 0.15, "nudge" = 0.225),
                           pal = NULL,
                           font_size = 11,
                           font = "",
                           ...) {

  if (is.null(pal)) {
    pal <- c("#ffc9b5", "#648767", "#b1ddf1", "#95a7ce", "#987284", "#3d5a80")
  }

  l <- list(...)
  if (is.null(l$alpha_par_nms)) alpha_par_nms <- NULL # just to show it's here
  if (is.null(l$rhat_upper)) rhat_upper <- 1.1
  if (is.null(l$ess_lower)) ess_lower <- 100
  if (is.null(l$moment_type)) moment_type <- "mean"

  type <- match.arg(type, c("parameter", "factor"))
  cred <- sort(cred)
  cred_l1 <- (1 - cred[2]) / 2
  cred_l2 <- (1 - cred[1]) / 2

  subjID <- value <- NULL # appease R CMD check

  all_data <- list()
  for (s in seq_along(summary_df)) {
    all_data[[s]] <- make_par_df(
      raw_df[[s]], summary_df[[s]], rhat_upper = rhat_upper,
      ess_lower = ess_lower
    )
  }
  all_data <- data.table::rbindlist(all_data, use.names = TRUE)
  if (is.null(pars_of_interest)) pars_of_interest <- unique(all_data[[par_nms]])
  all_data <- all_data[all_data[[par_nms]] %in% pars_of_interest, ]

  if (!is.null(afct_filt)) {
    all_data <- all_data[all_data[[afct_col]] == afct_filt, ]
  }

  if (!is.null(combine_with)) {
    all_data[[par_nms]] <- paste(
      all_data[[par_nms]], all_data[[combine_with]], sep = "_"
    )
  }
  if (!is.null(recode_na)) {
    all_data <- all_data |>
      dplyr::mutate(
        dplyr::across(.cols = by, .fns = ~ifelse(is.na(.), recode_na, .))
      )
  }

  if (!is.null(factor_scores)) {
    all_data <- all_data |>
      dplyr::left_join(factor_scores, by = "subjID")
  }
  if (is.null(par_lvls)) par_lvls <- unique(all_data[[par_nms]])
  if (flip) {
    par_lvls <- rev(par_lvls)
    if (!is.null(par_labs)) par_labs <- rev(par_labs)
  }

  if (is.null(par_labs)) par_labs <- par_lvls

  val_nm <- rlang::sym(par_vals)
  if (type == "parameter") {
    df <- all_data |> dplyr::rename(value = !!val_nm)
    df[[par_nms]] <- factor(df[[par_nms]], levels = par_lvls, labels = par_labs)
  } else if (type == "factor") {
    df <- all_data |>
      tidyr::pivot_longer(cols = par_lvls, names_to = type) |>
      dplyr::distinct(subjID, factor, .keep_all = TRUE)
    df[[type]] <- factor(df[[type]], levels = par_lvls, labels = par_labs)
  }

  type <- rlang::sym(type)

  if (is.null(by)) {
    rain_plot <- df |>
      ggplot2::ggplot(ggplot2::aes(x = !!type, y = value, fill = !!type,
                                   colour = !!type)) +
      ggplot2::guides(colour = "none", fill = "none") +
      ggplot2::scale_colour_manual(values = pal) +
      ggplot2::scale_fill_manual(values = pal)
  } else {
    by <- rlang::sym(by)
    by_length <- length(unique(df[[by]]))
    rain_plot <- df |>
      ggplot2::ggplot(ggplot2::aes(x = !!type, y = value, fill = factor(!!by),
                                   colour = interaction(!!by, !!type))) +
      ggplot2::guides(colour = "none") +
      ggplot2::scale_colour_manual(
        values = rep(pal[1:by_length], length(par_lvls))
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
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = as.numeric(!!type) - scatter[[2]]),
      position = ggplot2::position_jitter(width = scatter[[1]], height = 0),
      size = 0.6,
      alpha = 0.2
    ) +
    ggplot2::stat_summary(
      geom = "boxplot",
      fun.data = function(x) {
        setNames(
          quantile_hdi(
            x,
            c(cred_l1, cred_l2, 0.5, 1 - cred_l2, 1 - cred_l1),
            transform = FALSE, moment = moment_type
          ), c("ymin", "lower", "middle", "upper", "ymax")
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

  if (any(is.null(scatter))) rain_plot$layers <- rain_plot$layers[c(1, 3)]
  if (flip) rain_plot <- rain_plot + ggplot2::coord_flip()

  if (type == "parameter") {
    x_labels <- list()
    for (i in seq_along(par_lvls)) {
      par <- par_labs[[i]]
      alpha_par <- grepl("alpha", par)
      x_labels[[i]] <- bquote(
        .(rlang::parse_expr(axis_title(par, i, test, alpha_par, alpha_par_nms)))
      )
    }
    rain_plot <- rain_plot +
      ggplot2::scale_x_discrete(name = NULL, labels = x_labels)

    if ("beta" %in% par_lvls) {
      rain_plot <- rain_plot +
        ggplot2::scale_y_continuous(
          name = "Posterior mean", breaks = c(0, 0.5, 1, 2, 4, 8),
          trans = "pseudo_log"
        )
    } else {
      rain_plot <- rain_plot +
        ggplot2::scale_y_continuous(name = "Posterior mean")
    }
  } else if (type == "factor") {
    rain_plot <- rain_plot +
      ggplot2::scale_x_discrete(name = NULL, labels = par_labs) +
      ggplot2::ylab("Predicted factor score")
  }
  return(rain_plot)
}
