#' Plot an individual's mean posterior predictions compared to their raw affect
#' ratings
#'
#' \code{plot_affect} is capable of plotting either grouped or
#' individual-level posterior predictions (vs. raw observations) for a defined
#' list of posterior predictions and/or grouping.
#'
#' @param fit_list List of outputs from [get_affect_ppc].
#' @param plt_type Possible types are "grouped" or "individual".
#' @param adj_order Same as [fit_learning_model()].
#' @param nouns Formatted noun versions of the adjectives, in order.
#' @param id_no If \code{grouped == FALSE}, a participant number to plot. If
#' left as \code{NULL}, defaults to the individual with the median \eqn{R^2} for
#' each adjective.
#' @param r2_coords If \code{grouped == FALSE}, coordinates to print the
#' \eqn{R^2} value.
#' @param legend_pos,pal,font,font_size Same as [plot_import].
#'
#' @return A single or list of \code{ggplot} object(s) depending on type.
#'
#' @importFrom magrittr %>%
#' @importFrom stats quantile
#'
#' @examples \dontrun{
#' fit_affect <- fit_learning_model(
#'   example_data$nd,
#'   model = "2a",
#'   affect = TRUE,
#'   exp_part = "training"
#' )
#'
#' fit_dfs <- list()
#' for (adj in c("happy", "confident", "engaged")) {
#'   fits_dfs[[adj]] <- get_affect_ppc(fit_affect$draws, fit_affect$raw_df, adj = adj)
#' }
#'
#' # Grouped plot
#' plot_affect(fit_dfs, plt_type = "grouped")
#'
#' # Individual-level median posterior predictions
#' plot_affect(fit_dfs, plt_type = "individual", r2_coords = c(0.8, 0.97))
#' }
#'
#' @export

plot_affect <- function(fit_list,
                        plt_type = c("individual", "grouped"),
                        adj_order = c("happy", "confident", "engaged"),
                        nouns = c("Happiness", "Confidence", "Engagement"),
                        id_no = NULL,
                        r2_coords = c(0.9, 0.8),
                        legend_pos = "right",
                        pal = NULL,
                        font = "",
                        font_size = 11) {

  plt_type <- match.arg(plt_type)

  type <- trial_no_q <- value <- mean_val <- se_val <- se_pred <- NULL

  if (plt_type == "grouped") {
    if(is.null(pal)) pal <- c("#ffc9b5", "#95a7ce", "#987284")
    ppc_list <- lapply(
      1:length(adj_order),
      function(f) {
        dplyr::mutate(
          data.table::rbindlist(
            fit_list[[f]]$indiv_ppcs, idcol = "subjID"
          ),
          adj = adj_order[f]
        )
      }
    )

    grouped_plot <-
      data.table::rbindlist(ppc_list) %>%
      dplyr::group_by(adj, type, trial_no_q) %>%
      dplyr::mutate(mean_val = mean(value), se_val = std(value)) %>%
      dplyr::distinct(trial_no_q, adj, type, mean_val, se_val) %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = trial_no_q, y = mean_val, color = adj, fill = adj,
          linetype = factor(type, levels = c("raw", "pred")))
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0, 120), breaks = seq(0, 120, 20)) +
      ggplot2::scale_y_continuous(limits = c(25, 72)) +
      ggplot2::geom_line(size = 1.1, alpha = 0.5) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = mean_val - se_val, ymax = mean_val + se_val),
        alpha = 0.3, colour = NA
      ) +
      ggplot2::scale_color_manual(
        name = NULL,
        values = pal
      ) +
      ggplot2::scale_fill_manual(
        name = NULL,
        values = pal
      ) +
      ggplot2::xlab("Trial number") +
      ggplot2::ylab(paste0("Mean (\u00B1 SE) affect rating")) +
      ggplot2::guides(linetype = "none", color = "none", fill = "none") +
      cowplot::theme_half_open(
        font_size = font_size,
        font = font
      ) +
      ggplot2::theme(legend.position = c(0.85, 0.85))
    return(grouped_plot)
  }

  else if (plt_type == "individual") {
    if(is.null(pal)) {
      pal <- c("#ffc9b5", "#648767", "#b1ddf1", "#95a7ce", "#987284", "#3d5a80")
    }
    median_id <- function(df, kind, id = NULL) {
      if (kind == "num") {
        subset(df, R2 == quantile(df$R2, p = 0.5, type = 1, na.rm = T))$id_no
      } else if (kind == "id" & !is.null(id)) {
        subset(df, id_no == id)$subjID
      }
    }

    len <- length(fit_list)
    id_vec <- vector(mode = "integer", length = len)

    if (is.null(id_no)) {
      id_vec <- sapply(1:len,
                       function(f) median_id(fit_list[[f]]$fit_df, "num"))
    } else {
      id_vec <- rep(id_no, len)
    }

    indiv_ppc_plots <- list()

    for (a in 1:length(adj_order)) {
      adj <- adj_order[a]
      r2 <- subset(fit_list[[a]]$fit_df, id_no == id_vec[a])$R2

      indiv_ppc_plots[[adj]] <-
        fit_list[[a]]$indiv_ppcs[[
          median_id(fit_list[[a]]$fit_df, "id", id_vec[a])
        ]] %>%
        ggplot2::ggplot(
          ggplot2::aes(x = trial_no_q, y = value, color = type, fill = type)
        ) +
        ggplot2::geom_line() +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = value - se_pred, ymax = value + se_pred),
          alpha = 0.5
        ) +
        ggplot2::scale_color_manual(
          name = "Data type", labels = c("Predicted", "Real"),
          values = c(pal[a*2-1], pal[a*2])
        ) +
        ggplot2::scale_fill_manual(
          name = "Data type", labels = c("Predicted", "Real"),
          values = c(pal[a*2-1], pal[a*2])
        ) +
        ggplot2::xlab("Trial number") +
        ggplot2::ylab(paste0(nouns[a], " rating /100")) +
        cowplot::theme_half_open(font_size = font_size, font = font) +
        ggplot2::annotation_custom(
          grid::textGrob(
            bquote(R^2~"="~.(round(r2, 3))),
            gp = grid::gpar(fontsize = 16, col = "steelblue4"),
            x = r2_coords[1], y = r2_coords[2]
          )
        )  +
        ggplot2::theme(legend.position = legend_pos)
    }
    return(indiv_ppc_plots)
  }
}
