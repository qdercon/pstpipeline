#' Plot an individual's mean posterior predictions compared to their raw affect
#' ratings
#'
#' \code{plot_affect} is capable of plotting either grouped or
#' individual-level posterior predictions (vs. raw observations) for a defined
#' list of posterior predictions and/or grouping.
#'
#' @param data Either a list of outputs from [get_affect_ppc], or weights from
#' [get_affect_wts].
#' @param plt_type Possible types are "grouped" or "individual" (for
#' [get_affect_ppc] outputs) or "weights" (for [get_affect_wts] output).
#' @param adj_order Same as [fit_learning_model()].
#' @param nouns Formatted noun versions of the adjectives, in order.
#' @param id_no If \code{grouped == FALSE}, a participant number to plot. If
#' left as \code{NULL}, defaults to the individual with the median \eqn{R^2} for
#' each adjective.
#' @param r2_coords If \code{grouped == FALSE}, coordinates to print the
#' \eqn{R^2} value.
#' @param cred Same as [plot_glm], ignored unless \code{plt_type == "weights"}.
#' @param legend_pos,pal,font,font_size Same as [plot_import].
#'
#' @return A single or list of \code{ggplot} object(s) depending on type.
#'
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
#'
#' # Weight plot
#' aff_wts <- get_affect_wts(fit_affect$summary)
#' plot_affect(aff_wts, plt_type = "weights"))
#' }
#'
#' @export

plot_affect <- function(data,
                        plt_type = c("individual", "grouped", "weights"),
                        adj_order = c("happy", "confident", "engaged"),
                        nouns = c("Happiness", "Confidence", "Engagement"),
                        id_no = NULL,
                        r2_coords = c(0.9, 0.8),
                        cred = c(0.95, 0.99),
                        legend_pos = "right",
                        pal = NULL,
                        font = "",
                        font_size = 11) {

  plt_type <- match.arg(plt_type)

  type <- trial_no_q <- value <- mean_val <- se_val <- se_pred <- name <-
    param <- post_mean <- NULL

  if (plt_type == "weights") {

    cred <- sort(cred)
    cred_l1 <- (1-cred[2])/2
    cred_l2 <- (1-cred[1])/2

    p <- unique(data$param)
    labs <- c("\u03B3", expression(w[0]), expression(w[1]), expression(w[1]^o),
              expression(w[1]^b), expression(w[2]), expression(w[3]))

    if (any(grepl("w1_o", p) & !grepl("w1_b", p))) labs <- labs[c(1,2,3,6,7)]
    else if (!grepl("w1_o", p)) labs <- labs[c(1,2,6,7)]

    weight_plot <- data |>
      dplyr::filter(!grepl("mu|sigma", param)) |>
      dplyr::mutate(
        adj = paste0(toupper(substr(adj, 1, 1)), substr(adj, 2, nchar(name)))
      ) |>
      ggplot2::ggplot(
        ggplot2::aes(
          x = param, y = post_mean, color = factor(adj), fill = factor(adj)
        )) +
      geom_flat_violin(position = ggplot2::position_nudge(x = .125, y = 0),
                       adjust = 2, trim = FALSE, alpha = 0.5) +
      ggplot2::geom_point(
        ggplot2::aes(x = as.numeric(as.factor(param)) - 0.225),
        position = ggplot2::position_jitter(width = .1, height = 0),
        size = .25,
        alpha = 0.15
      ) +
      ggplot2::stat_summary(
        geom = "boxplot",
        fun.data = function(x) {
          setNames(
            quantile_hdi(x, c(cred_l1, cred_l2, 0.5, 1 - cred_l2, 1 - cred_l1)),
            c("ymin", "lower", "middle", "upper", "ymax")
          )
        },
        position = ggplot2::position_dodge2(),
        alpha = 0.6,
        width = 0.2
      ) +
      ggplot2::geom_hline(
        ggplot2::aes(yintercept = 0),
        size = 0.6,
        alpha = 0.4,
        linetype = "dashed",
        colour = "slategrey"
      ) +
      ggplot2::scale_color_manual(name = NULL, values = pal) +
      ggplot2::scale_fill_manual(name = NULL, values = pal) +
      ggplot2::scale_x_discrete(name = "Parameter", labels = labs) +
      ggplot2::scale_y_continuous(name = "Posterior mean") +
      cowplot::theme_half_open(font_family = font, font_size = font_size) +
      ggplot2::theme(legend.position = legend_pos)

  }
  else if (plt_type == "grouped") {
    if(is.null(pal)) pal <- c("#ffc9b5", "#95a7ce", "#987284")
    ppc_list <- lapply(
      1:length(adj_order),
      function(f) {
        dplyr::mutate(
          data.table::rbindlist(
            data[[f]]$indiv_ppcs, idcol = "subjID"
          ),
          adj = adj_order[f]
        )
      }
    )

    grouped_plot <-
      data.table::rbindlist(ppc_list) |>
      dplyr::group_by(adj, type, trial_no_q) |>
      dplyr::mutate(mean_val = mean(value), se_val = std(value)) |>
      dplyr::distinct(trial_no_q, adj, type, mean_val, se_val) |>
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

    len <- length(data)
    id_vec <- vector(mode = "integer", length = len)

    if (is.null(id_no)) {
      id_vec <- sapply(1:len,
                       function(f) median_id(data[[f]]$fit_df, "num"))
    } else {
      id_vec <- rep(id_no, len)
    }

    indiv_ppc_plots <- list()

    for (a in 1:length(adj_order)) {
      adj <- adj_order[a]
      r2 <- subset(data[[a]]$fit_df, id_no == id_vec[a])$R2

      indiv_ppc_plots[[adj]] <-
        data[[a]]$indiv_ppcs[[
          median_id(data[[a]]$fit_df, "id", id_vec[a])
        ]] |>
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
