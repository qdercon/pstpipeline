#' Plot associations between Q-learning model parameters and variables of interest
#'
#' \code{plot_glm} is designed to plot the results of fitting Bayesian GLMs, with
#' QL model parameters as the response variables, via the [parameter_glm()]
#' function. It plots both a box-and-whisker plot (defaulting to 95% HDIs for the box,
#' and 99% HDIs for the lines), plus the posterior distribution of coefficients, using
#' a version of [ggbeeswarm::geom_quasirandom()] which has been modified to output
#' half-violin plots made up of (half of) the individual posterior draws.
#'
#' @param par_df A [posterior::draws_df()], likely obtained from running
#' [parameter_glm()].
#' @param plot_var The variable of interest to plot (e.g., distanced vs non-distanced).
#' @param id.col The column that contains the QL model parameter names.
#' @param grp Optional group to plot separately on each plot, which should be the interaction
#' variable specified in [parameter_glm()].
#' @param grp_labs Optional labels for the groups defined by \code{grp}. It is recommended to
#' first run the function with this kept as \code{NULL} to make sure you label the correct
#' densities.
#' @param ovrll_title Title to set for the whole plot.
#' @param cred Vector, length 2, which defines the % HDI covered by the boxplot boxes and
#' lines respectively.
#' @param top_right,coord_flip Booleans; if \code{TRUE} the densities will be on the top
#' or the right depending on whether \code{coord_flip} is \code{TRUE} or \code{FALSE}
#' respectively.
#' @param dist_nudge,max_dist_width,point_alpha,point_size Control the position and size of
#' the density, or the transparency and size of the points that make it up respectively.
#' @param box_alpha,box_width,box_nudge Control the transparency, size, and position of the
#' summary boxplot.
#' @param pal,font_size,font Same as [plot_import()].
#'
#'
#' @importFrom stats setNames

plot_glm <- function(par_df, plot_var, id.col = "parameter", grp = id.col,
                     grp_labs = NULL, ovrll_title = NULL, cred = c(0.95, 0.99),
                     top_right = TRUE, coord_flip = TRUE, dist_nudge = 0,
                     max_dist_width = 0.5, point_alpha = 0.25, point_size = 0.75,
                     box_alpha = 0.6, box_width = 0.125, box_nudge = 0.1,
                     pal = NULL, font_size = 11, font = "") {

  if (is.null(pal)) {
    pal <- c("#ffc9b5", "#648767", "#b1ddf1", "#95a7ce", "#987284", "#3d5a80")
  }
  if (font != "") {
    extrafont::loadfonts(device = "win", quiet = TRUE)
  }

  cred <- sort(cred)
  cred_l1 <- (1-cred[2])/2
  cred_l2 <- (1-cred[1])/2

  parameter <- value <- NULL

  par_df <- par_df %>%
    dplyr::rename(value = plot_var)

  plots <- list()
  pars <- unique(par_df[[id.col]])
  if (y_grp != id.col & is.null(y_labs)) {
    y_labs <- levels(factor(par_df[[y_grp]]))
  }
  y_grp <- rlang::sym(y_grp)

  for (p in seq_along(pars)) {
    par <- pars[[p]]
    alpha <- grepl("_", par)
    voi <- rlang::sym(plot_var)

    par_df_tr <- par_df %>% dplyr::filter(parameter == par)
    if (alpha) {
      title <- "Estimated mean % difference in"
      plot <- par_df_tr %>%
        ggplot2::ggplot(
          ggplot2::aes(x = !!y_grp, y = (exp(value) - 1)*100, fill = id.col, colour = id.col)
          )
    }
    else {
      title <- "Estimated mean difference in"
      plot <- par_df_tr %>%
        ggplot2::ggplot(
          ggplot2::aes(x = !!y_grp, y = value, fill = id.col, colour = id.col)
        )
    }
    if (p==1 | !coord_flip) y_labels <- y_labs
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
      ggplot2::geom_hline(ggplot2::aes(yintercept = 0), alpha = 0.5, linetype = "dashed") +
      ggplot2::guides(colour = "none", fill = "none") +
      ggplot2::scale_colour_manual(values = pal[p]) +
      ggplot2::scale_fill_manual(values = pal[p]) +
      ggplot2::scale_x_discrete(name = NULL, labels = y_labels) +
      ggplot2::scale_y_continuous(
        name = bquote(
          .(title) ~ .(
            rlang::parse_expr(
              paste0(strsplit(par, "_")[[1]][1],
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
    title <- cowplot::ggdraw() +
      cowplot::draw_label(
        ovrll_title, x = 0, hjust = 0,
        fontfamily = font, size = 16,
        fontface = "bold"
      ) +
      ggplot2::theme(
        plot.margin = ggplot2::margin(0, 0, 0, 7)
      )
    plot_together <- cowplot::plot_grid(
      title, plot_together,
      nrow = 2, rel_heights = c(0.1, 1)
    )
  }
  return(plot_together)
}
