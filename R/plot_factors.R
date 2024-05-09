#' Various plots for factor prediction and derivation
#'
#' \code{plot_factors} plots various aspects of the transdiagnostic factor
#' derivation, including histograms and heatmaps of the questions themselves.
#'
#' @param df Data frame with factor scores or questions to plot. For
#' \code{predictive} and \code{factor_htmp} plots, this should be a named list
#' of \code{data.frames} (including "preds" & "scores" or "qns" & "coefs"
#' respectively). If a \code{factor_hist} plot with \code{grouped = TRUE} is
#' desired, it should be a list of \code{data.frames} (one per group).
#' @param plot_type Plot(s) to output: \code{factor_hist} (which can be
#' grouped), \code{r2_plot},\code{predictive}, and \code{factor_htmp}.
#' @param colnames Column names to get data from.
#' @param titles Title(s) for the plot(s) or pretty names more generally.
#' @param r2 Array of \eqn{R^{2}}{R^2} values.
#' @param qn Number of questions - used to find correct \eqn{R^{2}}{R^2} values.
#' @param hyp_alph Chosen alpha value (used to draw a dotted line on an
#' \code{r2_plot}).
#' @param grouped .
#' @param pal,font,font_size Same as [plot_import].
#'
#' @return A single or \code{list} of \code{ggplot} object(s).
#'
#' @examples
#' # See the notebook (data_cleaning_factor_derivation.ipynb) for examples as
#' # it requires Python input.
#' @export

plot_factors <- function(df,
                         plot_type,
                         colnames = NA,
                         titles = NA,
                         r2 = NA,
                         qn = NA,
                         hyp_alph = 0.1,
                         grouped = FALSE,
                         pal = NULL,
                         font = "",
                         font_size = 11) {
  if (is.null(pal)) {
    pal <- c("#ffc9b5", "#648767", "#b1ddf1", "#95a7ce", "#987284", "#3d5a80")
  }
  ret <- list()

  ## to appease R CMD check
  Factor <- Score <- Weight <- ..count.. <- alpha <- n_items <- id <- value <-
    obs <- predicted <- question <- NULL

  if (any(plot_type == "factor_hist")) {

    if (grouped) {
      all_datasets <- list()
      for (d in seq_along(df)) {
        all_datasets[[d]] <- df[[d]] |>
          dplyr::mutate(dataset = paste0("group_", d))
      }
      df <- data.table::rbindlist(all_datasets, use.names = TRUE)
      pal <- split(
        pal, ceiling(seq_along(pal) / length(unique(df[["dataset"]])))
      )
      group <- rlang::sym("dataset")
    } else {
      group <- rlang::sym("Factor")
    }
    hist_factors <- list()
    for (f in seq_along(colnames)) {
      hist_plot <- df |>
        tidyr::pivot_longer(
          cols = tidyselect::all_of(colnames), names_to = "Factor",
          values_to = "Score"
        ) |>
        dplyr::filter(Factor == colnames[f]) |>
        ggplot2::ggplot(ggplot2::aes(x = Score)) +
        ggplot2::geom_histogram(
          ggplot2::aes(y = ..count.., fill = !!group),
          colour = "black", alpha = 0.65, binwidth = 0.2, position = "identity"
        )  +
        ggplot2::geom_line(
          ggplot2::aes(y = ..count.., colour = !!group),
          binwidth = 0.2, stat = "bin"
        ) +
        ggplot2::scale_colour_manual(values = unlist(pal[[f]])) +
        ggplot2::scale_fill_manual(values = unlist(pal[[f]])) +
        ggplot2::guides(colour = "none", fill = "none") +
        ggplot2::scale_y_continuous(name = "Count") +
        cowplot::theme_half_open(
          font_size = font_size,
          font_family = font
        ) +
        ggplot2::ggtitle(titles[f])

      hist_factors[[f]] <- hist_plot
    }
    ret$hist_factors_all <- cowplot::plot_grid(
      plotlist = hist_factors, nrow = 1
    )
  }
  if (any(plot_type == "r2_plot")) {
    df <- as.data.frame(df) |>
      tibble::rownames_to_column(var = "Factor") |>
      tidyr::pivot_longer(-Factor, names_to = c("alpha", "n_items"),
                          names_sep = "_", values_to = "R2") |>
      dplyr::mutate(alpha = as.numeric(alpha)) |>
      dplyr::mutate(n_items = as.numeric(n_items))

    n_item_vec <- unique(sort(df$n_items, decreasing = FALSE))
    alphas <- unique(sort(df$alpha, decreasing = TRUE))
    index_alph <- which(alphas == hyp_alph)

    r2_plot <- df |>
      ggplot2::ggplot(ggplot2::aes(x = n_items, y = R2, colour = Factor)) +
      ggplot2::geom_line(size = 0.8) +
      ggplot2::scale_colour_manual(values = unlist(pal)) +
      ggplot2::scale_x_continuous(
        name = "No. questions", breaks = n_item_vec,
        sec.axis = ggplot2::dup_axis(
          labels = alphas, name = expression(alpha)
        )
      ) +
      ggplot2::scale_y_continuous(name = expression(R^2)) +
      ggplot2::geom_vline(
        ggplot2::aes(xintercept = n_item_vec[index_alph]),
        linetype = "dashed"
      ) +
      cowplot::theme_half_open(font_size = font_size, font_family = font) +
      ggplot2::theme(
        legend.position = c(0.75, 0.25),
        axis.text.x.top = ggplot2::element_text(angle = 90)
      )

    ret$r2_plot <- r2_plot
  }
  if (any(plot_type == "predictive")) {
    if (!is.list(df)) stop("Need 'preds' and 'scores' in a named list.")
    r2_col <- round(r2[, grep(qn, colnames(r2))], 3)
    pred_plots <- list()
    scores <- tibble::as_tibble(df[["scores"]]) |>
      dplyr::mutate(value = "obs")
    df_all <-
      dplyr::bind_cols(
        scores$id, as.data.frame(df[["preds"]]),
        .name_repair = ~make.names(c("id", colnames), unique = TRUE)
      ) |>
      dplyr::mutate(value = "predicted") |>
      dplyr::bind_rows(scores) |>
      tidyr::pivot_longer(cols = c(-id, -value),
                          values_to = "Score", names_to = "Factor") |>
      tidyr::pivot_wider(names_from = value, values_from = Score)

    for (f in seq_along(colnames)) {
      pred_plots[[f]] <- df_all |>
        dplyr::filter(Factor == colnames[f]) |>
        ggplot2::ggplot(ggplot2::aes(x = obs, y = predicted)) +
        ggplot2::geom_point(size = 2, alpha = 0.5, fill = pal[[f]],
                            colour = pal[[f]]) +
        ggplot2::geom_smooth(method = "lm", formula = "y~x", se = FALSE,
                             fill = pal[[f]], colour = pal[[f]]) +
        ggplot2::guides(colour = "none", fill = "none") +
        cowplot::theme_half_open(font_size = font_size, font_family = font) +
        ggplot2::xlab("True score") +
        ggplot2::ylab("Predicted score") +
        ggplot2::ggtitle(
          titles[f], subtitle = bquote(R^2 ~ "=" ~ .(r2_col[[f]]))
        )
    }
    ret$pred_plot <- cowplot::plot_grid(plotlist = pred_plots, nrow = 1)
  }
  if (any(plot_type == "factor_htmp")) {
    if (!is.list(df)) stop("Need 'qns' and 'coefs' in a named list.")
    names <- names(df[["qns"]][-1])
    heatmap <-
      tibble::as_tibble(
        df[["coefs"]], .name_repair = ~make.names(colnames, unique = TRUE)
      ) |>
      dplyr::mutate(question = names) |>
      dplyr::filter(!dplyr::if_all(.cols = 1:3, ~ . == 0)) |>
      tidyr::pivot_longer(cols = 1:3, names_to = "Factor",
                          values_to = "Weight") |>
      ggplot2::ggplot(ggplot2::aes(x = factor(question, levels = names),
                                   y = factor(Factor, levels = rev(colnames)),
                                   fill = Weight)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_distiller(palette = "Blues", direction = 1) +
      ggplot2::scale_y_discrete(name = NULL, labels = rev(titles)) +
      cowplot::theme_half_open(font_size = font_size, font_family = font) +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)
      )
    ret$heatmap <- heatmap
  }

  if (length(ret) == 1) return(ret[[1]])
  else return(ret)
}
