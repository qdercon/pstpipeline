#' Plot raw experiment data (single result)
#'
#' \code{plot_factors} plots various aspects of the transdiagnostic factor derivation, including
#' histograms and heatmaps of the questions themselves.
#'
#' @param df Data frame with factor scores or questions to plot.
#' @param plot_type Plot(s) to output: "factor_hist" (which can be grouped), .
#' @param colnames Column names to get data from.
#' @param titles Title(s) for the plot(s).
#' @param pal Custom colour palette to use.
#' @param group Grouping for each of the plots (if relevant).
#' @param font Use a custom font for the plots?
#' @param font_size Base plot font size.
#'
#' @return \code{list} of \code{ggplot} objects
#'
#' @importFrom magrittr %>%
#' @export

plot_factors <- function(df, plot_type, colnames, titles, pal = NULL, group = NULL,
                         font_size = 11, font = NULL) {
  if (!is.null(font)) {
    extrafont::loadfonts(device = "win", quiet = TRUE)
  }
  if (is.null(pal)) {
    message("No colour palette selected (pal=NULL), reverting to defaults.")
    pal <- c('#ffa630', '#42bfdd', '#ef3e36', '#745296','#f08080','#d17a22')
  }

  ret <- list()

  if (any(plot_type == "factor_hist")) {
    hist_factors <- list()
    if (!is.null(group)) {
      pal <- split(pal, ceiling(seq_along(pal)/length(unique(df[[group]]))))
      group <- rlang::sym(group)
    }
    for (f in seq_along(colnames)) {
      hist_plot <- df %>%
        tidyr::pivot_longer(cols = colnames, names_to = "Factor", values_to = "Score") %>%
        dplyr::filter(Factor == colnames[f]) %>%
        ggplot2::ggplot(ggplot2::aes(x = Score))

      if (is.null(group)) {
        hist_plot <- hist_plot +
          ggplot2::geom_histogram(
            ggplot2::aes(y = ..count.., fill = Factor), colour = "black", alpha = 0.4,
                         binwidth = 0.2, position = "identity"
            ) +
          ggplot2::geom_line(
            ggplot2::aes(y = (..density..*(dim(factor_score_pred)[1]*0.2))),
            size = 1, stat = 'density', colour = pal[[f]])
      } else {
        hist_plot <- hist_plot +
          ggplot2::geom_histogram(
            ggplot2::aes(y = ..count.., fill = !!group), colour = "black", alpha = 0.4,
                         binwidth = 0.2, position = "identity"
            ) +
          ggplot2::geom_line(
            ggplot2::aes(y = (..density..*(dim(factor_score_pred)[1]*0.2)), colour = !!group),
            size = 1, stat = 'density'
            ) +
         ggplot2::scale_colour_manual(values = unlist(pal[[f]])) +
         ggplot2::guides(colour = "none")
      }

      hist_plot <- hist_plot +
        ggplot2::scale_fill_manual(values = pal[[f]]) +
        ggplot2::guides(fill = "none") +
        ggplot2::scale_y_continuous(name = "Count") +
        cowplot::theme_half_open(
          font_size = font_size,
          font_family = ifelse(!is.null(font), font, "")
        ) +
        ggplot2::ggtitle(titles[f])

      hist_factors[[f]] <- hist_plot
    }
    ret$hist_factors_all <- ggpubr::ggarrange(plotlist = hist_factors, nrow = 1)
  }

  if (length(ret)==1) return(ret[[1]])
  else return(ret)
}
