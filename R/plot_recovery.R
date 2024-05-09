#' Plots to check parameter recovery
#'
#' \code{plot_recovery} produces correlation plots between input and output
#' parameter values, plus a confusability matrix.
#'
#' @param raw_pars Parameters outputted from [simulate_QL()].
#' @param rec_pars [cmdstanr::summary()] containing the fitted parameter values
#' for the relevant model.
#' @param test Boolean indicating whether recovered parameters are from the test
#' phase.
#' @param affect Boolean indicating whether recovered parameters are from an
#' affect model.
#' @param alpha_par_nms Option to rename learning rate parameters, defaults to
#' the names from \code{par_df}.
#' @param plot_together If \code{TRUE}, returns a panel with all plots plotted,
#' otherwise a named list of plots is returned.
#' @param incl_legend Boolean indicating whether to include legends for affect
#' parameter plots (for the question types).
#' @param pal,font,font_size Same as [plot_import()].
#'
#' @return A named \code{list} of \code{ggplot} objects.
#'
#' @importFrom rlang !!
#'
#' @returns A list or grid of ggplots.
#'
#' @examples \dontrun{
#' dir.create("outputs/cmdstan/simulated_data")
#'
#' train_sim_2a <- simulate_QL(
#'   sample_size = 10,
#'   alpha_pos_dens = c(shape = 2, scale = 0.1), # default
#'   alpha_neg_dens = c(shape = 2, scale = 0.1), # default
#'   beta_dens = c(mean = 3, sd = 1) # default
#' )
#'
#' mcmc_2a_train_sim <- fit_learning_model(
#'   train_sim_2a$sim, model = "2a", exp_part = "training",
#'   vb = FALSE, model_checks = FALSE, par_recovery = TRUE,
#'   outputs = c("raw_df", "stan_datalist", "summary", "draws_list"),
#'   out_dir = "outputs/cmdstan/simulated_data"
#' )
#'
#' plot_recovery(train_sim_2a$pars, mcmc_2a_train_sim$summary)
#' }
#'
#' @export

plot_recovery <- function(raw_pars,
                          rec_pars,
                          test = FALSE,
                          affect = FALSE,
                          alpha_par_nms = NULL,
                          plot_together = TRUE,
                          incl_legend = !plot_together,
                          pal = NULL,
                          font = "",
                          font_size = 11) {

  if (is.null(pal)) {
    pal <- c("#ffc9b5", "#648767", "#b1ddf1", "#95a7ce", "#987284", "#3d5a80")
  }

  # to appease R CMD check
  parameter <- obs_mean <- rec_mean <- id_no <- rec_var <- obs_var <- corr <-
    aff_num <- adj <- outc_lag <- NULL

  rec_pars_df <- clean_summary(rec_pars) |>
    dplyr::rename(rec_mean = mean)

  pars_df <- raw_pars |>
    dplyr::select(-tidyselect::any_of("subjID")) |>
    tidyr::pivot_longer(
      cols = c(tidyselect::matches("alpha|beta|gamma|w")),
      names_to = "parameter", values_to = "obs_mean"
    ) |>
    dplyr::left_join(rec_pars_df) |>
    suppressMessages() |>
    dplyr::select(-tidyselect::any_of(c("rhat", "ess_bulk", "ess_tail")))

  pars <- list()

  if (affect) {
    pars_df <- pars_df |> tidyr::drop_na(obs_mean, rec_mean)
    if ("outc_lag" %in% colnames(pars_df)) {
      pars_df <- pars_df |>
        dplyr::mutate(
          parameter = ifelse(
            !is.na(outc_lag), paste(parameter, outc_lag, sep = "_"), parameter
          )
        ) |>
        #dplyr::arrange(parameter) |>
        dplyr::select(-outc_lag)
    }
    pars[[1]] <- pars_df |> dplyr::filter(is.na(aff_num) & !is.na(obs_mean))
    pars[[2]] <- pars_df |>
      dplyr::filter(!is.na(aff_num) & !is.na(obs_mean)) |>
      dplyr::mutate(
        adj = factor(paste0(
          toupper(substr(adj, 1, 1)), substr(adj, 2, nchar(adj))
        ))
      )
  } else {
    pars[[1]] <- pars_df
  }

  pred_plots <- list()
  pred_plots$cor_plots <- list()
  all_pars <- unique(pars_df$parameter)
  num_ql <- sum(grepl("alpha|beta", all_pars))

  if (length(all_pars) > length(pal)) {
    pal <- rep(pal, 1 + length(all_pars) %/% length(pal))
  }

  for (p in seq_along(all_pars)) {
    legend_pos <- "none"
    if (!affect || p <= num_ql) {
      pars_to_plot <- pars[[1]] |> dplyr::mutate(adj = factor(parameter))
      line_col <- col_vals <- pal[[p]]
    } else {
      pars_to_plot <- pars[[2]]
      line_col <- "slategray"
      col_vals <- pal[(num_ql + 1):length(pal)]
      if (incl_legend) legend_pos <- "right"
    }

    par <- all_pars[p]
    alpha_par <- grepl("alpha", par)

    pred_plots$cor_plots[[all_pars[p]]] <- pars_to_plot |>
      dplyr::filter(parameter == par) |>
      ggplot2::ggplot(
        ggplot2::aes(x = obs_mean, y = rec_mean, fill = adj, colour = adj)
      ) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_smooth(
        method = "lm", formula = "y~x", se = FALSE, fill = line_col,
        colour = line_col
      ) +
      ggplot2::scale_colour_manual(name = NULL, values = col_vals) +
      ggplot2::scale_fill_manual(name = NULL, values = col_vals) +
      ggplot2::ggtitle(
        bquote(
          .(rlang::parse_expr(
            axis_title(par, p, test, alpha_par, alpha_par_nms)
          ))
        ),
        subtitle = bquote(
          r ~ "=" ~ .(
            round(
              cor(pars_to_plot[pars_to_plot$parameter == par, ]$obs_mean,
                pars_to_plot[pars_to_plot$parameter == par, ]$rec_mean
              ), 2
            )
          )
        )
      ) +
      ggplot2::xlab("Observed") +
      ggplot2::ylab("Recovered") +
      cowplot::theme_half_open(
        font_size = font_size,
        font_family = font
      ) +
      ggplot2::theme(legend.position = legend_pos)
  }

  htmps <- list()

  for (c in seq_along(pars)) {
    cor_mat <- pars[[c]] |>
      tidyr::pivot_wider(names_from = parameter,
                         values_from = tidyselect::contains("mean")) |>
      dplyr::select(
        -tidyselect::any_of(tidyselect::matches("aff_num|adj")),
      ) |>
      dplyr::select(-id_no) |>
      cor()

    par_nms <- unique(pars[[c]]$parameter)
    labs <- sapply(
      seq_along(par_nms),
      function(n) {
        axis_title(
          par_nms[n], n, test, grepl("alpha", par_nms[n]), alpha_par_nms
        )
      }
    )

    htmps[[c]] <-
      tibble::as_tibble(cor_mat, rownames = "rec_var") |>
      dplyr::select(tidyselect::matches("var|obs")) |>
      dplyr::filter(grepl("rec", rec_var)) |>
      dplyr::mutate(rec_var = sub("rec_mean_", "", rec_var)) |>
      tidyr::pivot_longer(cols = tidyselect::contains("obs"),
                          names_to = "obs_var", names_prefix = "obs_mean_",
                          values_to = "corr") |>
      dplyr::mutate(rec_var = factor(rec_var, levels = par_nms)) |>
      dplyr::mutate(obs_var = factor(obs_var, levels = par_nms)) |>
      ggplot2::ggplot(ggplot2::aes(x = obs_var, y = rec_var, fill = corr)) +
      ggplot2::geom_tile(alpha = 0.6) +
      ggplot2::guides(fill = "none") +
      ggplot2::geom_text(ggplot2::aes(label = round(corr, 2), family = font)) +
      ggplot2::scale_fill_distiller(palette = "Blues", direction = 1) +
      cowplot::theme_half_open(font_size = font_size, font_family = font) +
      ggplot2::scale_x_discrete(
        name = "Observed", labels = rlang::parse_exprs(labs)
      ) +
      ggplot2::scale_y_discrete(
        name = "Recovered", labels = rlang::parse_exprs(labs)
      )
  }

  if (!affect) {
    pred_plots$heatmap <- htmps[[1]]
    if (!plot_together) {
      return(pred_plots)
    } else {
      ql_plots <- pred_plots$cor_plots
      ql_plots$heatmap <- pred_plots$heatmap
      plot <- cowplot::plot_grid(plotlist = ql_plots, nrow = 1)
      return(plot)
    }
  } else {
    pred_plots$heatmaps <- list()
    pred_plots$heatmaps$ql <- htmps[[1]]
    pred_plots$heatmaps$wts <- htmps[[2]]
    if (!plot_together) {
      return(pred_plots)
    } else {
      ql_plots <- pred_plots$cor_plots[grep("alpha|beta", all_pars)]
      ql_plots$heatmap <- pred_plots$heatmap$ql
      wt_plots <- pred_plots$cor_plots[grep("w|gamma", all_pars)]
      wt_plots$heatmap <- pred_plots$heatmap$wts

      nwts <- length(grep("w", all_pars))

      if (nwts == 3) {
        plot <- cowplot::plot_grid(
          cowplot::plot_grid(plotlist = ql_plots, nrow = 1),
          cowplot::plot_grid(plotlist = wt_plots, nrow = 1),
          nrow = 2
        )
      } else if (nwts == 4) {
        plot <- cowplot::plot_grid(
          cowplot::plot_grid(plotlist = ql_plots, nrow = 1),
          cowplot::plot_grid(
            plotlist = wt_plots[c(1, 2, 6, 3, 4, 5)], nrow = 2
          ),
          nrow = 2,
          rel_heights = c(1, 2)
        )
      } else if (nwts == 5) {
        plot <- cowplot::plot_grid(
          cowplot::plot_grid(plotlist = ql_plots, nrow = 1),
          cowplot::plot_grid(
            cowplot::plot_grid(
              plotlist = wt_plots[-length(wt_plots)], nrow = 2
            ),
            wt_plots$heatmap,
            ncol = 2,
            rel_widths = c(3, 1)
          ),
          nrow = 2,
          rel_heights = c(1, 2)
        )
      } else {
        message("Unable to coerce plots together, returning list of plots.")
        return(pred_plots)
      }
      return(plot)
    }
  }
}
