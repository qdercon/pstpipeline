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
#' @param grp_labs Optional labels for the groups defined by \code{grp}. It is
#' recommended to first run the function with this kept as \code{NULL} to make
#' sure you label the correct densities.
#' @param ovrll_title Title to set for the whole plot.
#' @param title_font_size,title_rel_ht Font size, and relative height of the
#' title compared to the main plot, given as a 2 element vector. Defaults to
#' \code{16}pt and \code{c(0.15, 1)} respectively.
#' @param plt_rel_widths Relative widths for the plotted grid (passed to
#' [cowplot::plot_grid()]). Defaults to \code{1} (i.e. all equal).
#' @param cred Vector, length 2, which defines the % HDI covered by the boxplot
#' boxes and lines respectively.
#' @param top_right,coord_flip Booleans; if \code{TRUE} the densities will be
#' on the top or the right depending on whether \code{coord_flip} is \code{TRUE}
#' or \code{FALSE} respectively.
#' @param dist_nudge,max_dist_width,point_alpha,point_size Control the position
#' and size of the density, or the transparency and size of the points that make
#' it up respectively.
#' @param box_alpha,box_width,box_nudge Control the transparency, size, and
#' position of the summary boxplot.
#' @param pal,font_size,font Same as [plot_import()].
#'
#' @importFrom stats setNames
#' @export

plot_glm <- function(
  par_df,
  plot_var,
  id.col = "parameter",
  test = FALSE,
  grp = id.col,
  grp_labs = NULL,
  ovrll_title = NULL,
  title_font_size = 16,
  title_rel_ht = NULL,
  plt_rel_widths = 1,
  cred = c(0.95, 0.99),
  top_right = TRUE,
  coord_flip = TRUE,
  dist_nudge = 0,
  max_dist_width = 0.5,
  point_alpha = 0.25,
  point_size = 0.75,
  box_alpha = 0.6,
  box_width = 0.125,
  box_nudge = 0.1,
  pal = NULL,
  font_size = 11,
  font = ""
  )
{

  if (is.null(pal)) {
    pal <- c("#ffc9b5", "#648767", "#b1ddf1", "#95a7ce", "#987284", "#3d5a80")
  }

  cred <- sort(cred)
  cred_l1 <- (1-cred[2])/2
  cred_l2 <- (1-cred[1])/2

  parameter <- value <- NULL ## appease R CMD check

  par_df <- par_df %>%
    dplyr::rename(value = plot_var)

  plots <- list()
  pars <- unique(par_df[[id.col]])
  if (grp != id.col) {
    if (!is.null(par_df[[paste0(grp, "_recode")]])) {
      grp <- paste0(grp, "_recode")
    }
    if (is.null(grp_labs)) grp_labs <- levels(factor(par_df[[grp]]))
  }
  grp <- rlang::sym(grp)

  for (p in seq_along(pars)) {
    par <- pars[[p]]
    alpha <- grepl("alpha", par)
    voi <- rlang::sym(plot_var)

    par_df_tr <- par_df %>% dplyr::filter(parameter == par)
    if (alpha) {
      title <- "Estimated mean % difference in"
      plot <- par_df_tr %>%
        ggplot2::ggplot(
          ggplot2::aes(x = !!grp, y = (exp(value) - 1)*100, fill = id.col,
                       colour = id.col)
          )
    }
    else {
      title <- "Estimated mean difference in"
      plot <- par_df_tr %>%
        ggplot2::ggplot(
          ggplot2::aes(x = !!grp, y = value, fill = id.col, colour = id.col)
        )
    }
    if (p==1 | !coord_flip) y_labels <- grp_labs
    else y_labels <- NULL

    plots[[par]] <- plot +
      geom_quasirandom(
        alpha = point_alpha, size = point_size, width = max_dist_width,
        half = TRUE, right = top_right, nudge = dist_nudge
      ) +
      ggplot2::stat_summary(
        geom = "boxplot",
        fun.data = function(x) {
          setNames(
            quantile_hdi(
              x, c(cred_l1, cred_l2, 0.5, 1 - cred_l2, 1 - cred_l1),
              transform = alpha), c("ymin", "lower", "middle", "upper", "ymax")
            )
        },
        position = ggplot2::position_nudge(x = -box_nudge),
        alpha = box_alpha,
        width = box_width
      ) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = 0), alpha = 0.5,
                          linetype = "dashed") +
      ggplot2::guides(colour = "none", fill = "none") +
      ggplot2::scale_colour_manual(values = pal[p]) +
      ggplot2::scale_fill_manual(values = pal[p]) +
      ggplot2::scale_x_discrete(name = NULL, labels = y_labels) +
      ggplot2::scale_y_continuous(
        name = bquote(
          .(title) ~ .(
            rlang::parse_expr(
              paste0(
                strsplit(par, "_")[[1]][1], ifelse(test, "*minute", ""),
                ifelse(!alpha, "", paste0("[", strsplit(par, "_")[[1]][2], "]"))
              )
            )
          )
        )
      ) +
      cowplot::theme_half_open(
        font_size = font_size,
        font_family = font
      )

    if (coord_flip) plots[[par]] <- plots[[par]] + ggplot2::coord_flip()
  }

  plot_together <- cowplot::plot_grid(
    plotlist = plots,
    nrow = 1
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
    plot_together <- cowplot::plot_grid(
      title, plot_together,
      nrow = 2, rel_heights = title_rel_ht
    )
  }
  return(plot_together)
}
