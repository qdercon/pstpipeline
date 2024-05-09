#' Plot associations between Q-learning model parameters and variables of
#' interest
#'
#' \code{plot_glm} is designed to plot the results of fitting Bayesian GLMs,
#' with QL model parameters as the response variables, via the [parameter_glm()]
#' function. It plots both a box-and-whisker plot (defaulting to 95% HDIs for
#' the box, and 99% HDIs for the lines), plus the posterior distribution of
#' coefficients (half-violin plots made up of the individual posterior draws).
#'
#' @param par_df A [posterior::draws_df()], likely obtained from running
#' [parameter_glm()].
#' @param plot_var The variable of interest to plot (e.g., distanced vs
#' non-distanced).
#' @param id.col The column that contains the QL model parameter names.
#' @param test Boolean indicating whether summaries are from the test phase.
#' @param grp Optional group to plot separately on each plot, which should be
#' the interaction variable specified in [parameter_glm()].
#' @param fclr To what variable should the colour scheme be applied? Defaults to
#' \code{grp} = \code{id.col}.
#' @param axis_fixed Logical indicating whether the y-axis should be fixed
#' across all plots. Defaults to \code{FALSE}.
#' @param grp_labs Optional labels for the groups defined by \code{grp}. It is
#' recommended to first run the function with this kept as \code{NULL} to make
#' sure you label the correct densities.
#' @param plot_together If \code{TRUE}, returns a panel with all plots plotted
#' as defined by subsequent arguments. Otherwise a named list of plots is
#' returned.
#' @param ovrll_title Title to set for the whole plot. Ignored if
#' \code{plot_together == FALSE}.
#' @param title_font_size,title_rel_ht Font size, and relative height of the
#' title compared to the main plot, given as a 2 element vector. Defaults to
#' \code{16}pt and \code{c(0.15, 1)} respectively. Ignored if
#' \code{plot_together == FALSE}.
#' @param plt_rows,plt_rel_widths Number of rows, and relative widths for the
#' plotted grid (passed to [cowplot::plot_grid()]). Defaults to \code{1}
#' (i.e. single row, all equal widths). Ignored if \code{!plot_together}.
#' @param cred Vector, length 2, which defines the % HDI covered by the boxplot
#' boxes and lines respectively.
#' @param coord_flip Plot horizontal (\code{TRUE}) or vertical (\code{FALSE})
#' densities.
#' @param box_alpha,box_width,box_nudge Control the transparency, size, and
#' position of the summary boxplot.
#' @param pal,font_size,font Same as [plot_import()].
#' @param ... Additional arguments such as \code{alpha_par_nms} to alter axis
#' titles or control [quantile_hdi()].
#'
#' @importFrom stats setNames
#'
#' @returns Either a list of ggplots or a [cowplot::plot_grid()] if
#' \code{plot_together == TRUE}.
#'
#' @examples \dontrun{
#' Comparing parameters across groups
#'
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
#' distanced <- parameter_glm(
#'   summary_df = list(fit_nd$summary, fit_dis$summary),
#'   raw_df = list(fit_nd$raw_df, fit_dis$raw_df),
#'   var_of_interest = "distanced",
#'   covariates = c("age", "sex", "digit_span"),
#'   iter_warmup = 1000, iter_sampling = 1000
#' )
#' plot_glm(distanced, plot_var = "distanced")
#' }
#'
#' @export

plot_glm <- function(par_df,
                     plot_var,
                     id.col = "parameter",
                     test = FALSE,
                     grp = id.col,
                     fclr = id.col,
                     axis_fixed = FALSE,
                     grp_labs = NULL,
                     plot_together = TRUE,
                     ovrll_title = NULL,
                     title_font_size = 16,
                     title_rel_ht = NULL,
                     plt_rows = 1,
                     plt_rel_widths = 1,
                     cred = c(0.95, 0.99),
                     coord_flip = TRUE,
                     box_alpha = 0.6,
                     box_width = 0.125,
                     box_nudge = 0.1,
                     pal = NULL,
                     font_size = 11,
                     font = "",
                     ...) {

  if (is.null(pal)) {
    pal <- c("#ffc9b5", "#648767", "#b1ddf1", "#95a7ce", "#987284", "#3d5a80")
  }

  l <- list(...)
  if (is.null(l$alpha_par_nms)) l$alpha_par_nms <- NULL # just to show it's here
  if (is.null(l$moment_type)) moment_type <- "mean"

  if (!id.col %in% colnames(par_df)) {
    stop(paste0(id.col, " column not found in parameter data frame."))
  }

  cred <- sort(cred)
  cred_l1 <- (1 - cred[2]) / 2
  cred_l2 <- (1 - cred[1]) / 2

  value <- NULL ## appease R CMD check

  par_df <- par_df |>
    dplyr::rename(value = tidyselect::all_of(plot_var))

  plots <- list()
  pars <- unique(par_df[[id.col]])
  if (!is.factor(par_df[[grp]])) par_df[[grp]] <- factor(par_df[[grp]])
  if (grp != id.col) {
    if (!is.null(par_df[[paste0(grp, "_recode")]])) {
      grp <- paste0(grp, "_recode")
    }
    if (is.null(grp_labs)) grp_labs <- levels(factor(par_df[[grp]]))
  }
  grp <- rlang::sym(grp)
  nc <- length(unique(par_df[[fclr]]))
  if (fclr != id.col) pal <- lapply(seq_along(pars), function(x) rev(pal[1:nc]))
  fclr <- rlang::sym(fclr)
  id_col <- rlang::sym(id.col)

  if (axis_fixed) {
    min <- min(par_df[["value"]])
    max <- max(par_df[["value"]])
    ylims <- c(min * 1.05, max * 1.05)
  }

  for (p in seq_along(pars)) {
    par <- pars[[p]]
    plots[[par]] <- local({ # inelegant way of making sure plots don't overwrite
      par <- pars[[p]]
      alpha_par <- grepl("alpha", par)
      gamma_dist <- grepl("alpha|gamma", par)

      par_df_tr <- par_df |> dplyr::filter(!!id_col == par)
      if (gamma_dist) {
        title <- "Estimated mean % difference in"
        plot <- par_df_tr |>
          ggplot2::ggplot(
            ggplot2::aes(
              x = !!grp, y = (exp(value) - 1) * 100, fill = !!fclr,
              colour = !!fclr
            )
          )
      } else {
        title <- "Estimated mean difference in"
        plot <- par_df_tr |>
          ggplot2::ggplot(
            ggplot2::aes(x = !!grp, y = value, fill = !!fclr, colour = !!fclr)
          )
      }
      if (p == 1 | !coord_flip) y_labels <- grp_labs
      else y_labels <- NULL

      axs_ttl <- bquote(
        .(title) ~ .(
          rlang::parse_expr(
            axis_title(par, p, test, alpha_par, l$alpha_par_nms)
          )
        )
      )

      plot <- plot +
        geom_flat_violin() +
        ggplot2::stat_summary(
          geom = "boxplot",
          fun.data = function(x) {
            setNames(
              quantile_hdi(
                x, c(cred_l1, cred_l2, 0.5, 1 - cred_l2, 1 - cred_l1),
                transform = gamma_dist, moment = moment_type
              ),
              c("ymin", "lower", "middle", "upper", "ymax")
            )
          },
          position = ggplot2::position_nudge(x = -box_nudge),
          alpha = box_alpha,
          width = box_width
        ) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = 0), alpha = 0.5,
                            linetype = "dashed") +
        ggplot2::guides(colour = "none", fill = "none") +
        ggplot2::scale_colour_manual(values = pal[[p]]) +
        ggplot2::scale_fill_manual(values = pal[[p]]) +
        ggplot2::scale_x_discrete(name = NULL, labels = y_labels) +
        cowplot::theme_half_open(
          font_size = font_size,
          font_family = font
        )

      if (axis_fixed) {
        plot <- plot +
          ggplot2::scale_y_continuous(name = axs_ttl, limits = ylims)
      } else {
        plot <- plot +
          ggplot2::scale_y_continuous(name = axs_ttl)
      }

      if (coord_flip) plot <- plot + ggplot2::coord_flip()
      return(plot)
    })
  }

  if (!plot_together) {
    return(plots)
  } else {
    plot_panel <- cowplot::plot_grid(
      plotlist = plots, rel_widths = plt_rel_widths, nrow = plt_rows
    )

    if (!is.null(ovrll_title)) {
      if (is.null(title_rel_ht)) title_rel_ht <- c(0.15, 1)
      title <- cowplot::ggdraw() +
        cowplot::draw_label(
          ovrll_title, x = 0, hjust = 0,
          fontfamily = font, size = title_font_size,
          fontface = "bold"
        ) +
        ggplot2::theme(
          plot.margin = ggplot2::margin(2.5, 0, 2.5, 7)
        )
      plot_panel <- cowplot::plot_grid(
        title, plot_panel, nrow = 2, rel_heights = title_rel_ht
      )
    }
    return(plot_panel)
  }
}
