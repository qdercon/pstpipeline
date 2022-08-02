#' Plots to check parameter recovery
#'
#' \code{plot_recovery} produces correlation plots between input and output
#' parameter values, plus a confusability matrix.
#'
#' @param raw_pars Data frame with factor scores or questions to plot.
#' @param sim_pars [cmdstanr::summary()] containing the fitted parameter values
#' for the relevant model.
#' @param test Boolean indicating whether recovered parameters are from the test
#' phase.
#' @param alpha_par_nms Option to rename learning rate parameters, defaults to
#' the names from \code{par_df}.
#' @param pal,font,font_size Same as [plot_import()].
#'
#' @return A named \code{list} of \code{ggplot} objects.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @export

plot_recovery <- function(raw_pars,
                          sim_pars,
                          test = FALSE,
                          alpha_par_nms = NULL,
                          pal = NULL,
                          font = "",
                          font_size = 11) {

  if (is.null(pal)) {
    pal <- c("#ffc9b5", "#648767", "#b1ddf1", "#95a7ce", "#987284", "#3d5a80")
  }

  # to appease R CMD check
  variable <- . <- parameter <- obs_mean <- sim_mean <- id_no <- sim_var <-
    obs_var <- corr <- NULL

  sim_pars_df <- sim_pars %>%
    dplyr::filter(grepl("alpha|beta", variable)) %>%
    dplyr::filter(!grepl("_pr", variable)) %>%
    dplyr::filter(!grepl("mu_", variable)) %>%
    dplyr::select(variable, mean) %>%
    dplyr::mutate(
      id_no = as.numeric(sub("\\].*$", "", sub(".*\\[", "", .[["variable"]])))
      ) %>%
    dplyr::mutate(variable = sub("\\[.*$", "", .[["variable"]])) %>%
    dplyr::rename(parameter = variable) %>%
    dplyr::rename(sim_mean = mean)

  pars_df <- raw_pars %>%
    dplyr::select(-tidyselect::any_of("subjID")) %>%
    tidyr::pivot_longer(
      cols = c(tidyselect::contains("alpha"), tidyselect::contains("beta")),
      names_to = "parameter", values_to = "obs_mean"
      ) %>%
    dplyr::left_join(sim_pars_df, by = c("id_no", "parameter"))

  pred_plots <- list()

  for (p in seq_along(unique(pars_df$parameter))) {
    par <- unique(pars_df$parameter)[p]
    alpha <- grepl("alpha", par)

    pred_plots[[p]] <- pars_df %>%
      dplyr::filter(parameter == par) %>%
      ggplot2::ggplot(ggplot2::aes(x = obs_mean, y = sim_mean)) +
      ggplot2::geom_point(size= 2, alpha = 0.5, fill = pal[[p]],
                          colour = pal[[p]]) +
      ggplot2::geom_smooth(
        method = "lm", formula = "y~x", se = FALSE, fill = pal[[p]],
        colour = pal[[p]]
        ) +
      ggplot2::scale_colour_manual(values = pal[p]) +
      ggplot2::scale_fill_manual(values = pal[p]) +
      ggplot2::guides(colour = "none", fill = "none") +
      ggplot2::ggtitle(
        bquote(
          .(rlang::parse_expr(
            axis_title(par, p, test, alpha_par, alpha_par_nms)
            ))
        ),
        subtitle = bquote(
          r~"="~.(round(cor(pars_df[pars_df$parameter == par,]$obs_mean,
                              pars_df[pars_df$parameter == par,]$sim_mean), 2))
        )
      ) +
      ggplot2::xlab("Observed") +
      ggplot2::ylab("Recovered") +
      cowplot::theme_half_open(
        font_size = font_size,
        font_family = font
      )
  }

  cor_mat <- pars_df %>%
    tidyr::pivot_wider(names_from = parameter,
                       values_from = tidyselect::contains("mean")) %>%
    dplyr::select(-id_no) %>%
    cor()

  n_pars <- dim(cor_mat)[1]/2
  if (n_pars == 3) {
    order <- c("alpha_pos", "alpha_neg", "beta")
    if (is.null(alpha_par_nms)) alpha_par_nms <- c("pos", "neg")
    if (test) {
      labs <- c(
        paste0("alpha*minute[", alpha_par_nms[1], "]"),
        paste0("alpha*minute[", alpha_par_nms[2], "]"),
        "beta*minute")
    }
    else {
      labs <- c(
        paste0("alpha[", alpha_par_nms[1], "]"),
        paste0("alpha[", alpha_par_nms[2], "]"),
        "beta")
    }
  }
  else {
    order <- labs <- c("alpha", "beta")
    if (test) labs <- c("alpha*minute", "beta*minute")
  }

  corr_htmp <- tibble::as_tibble(
      cor_mat[-(1:n_pars), -((n_pars+1):(2*n_pars))], rownames = NA
    ) %>%
    dplyr::mutate(sim_var = sub("sim_mean_", "", row.names(.))) %>%
    tidyr::pivot_longer(cols = tidyselect::contains("obs"),
                        names_to = "obs_var", names_prefix = "obs_mean_",
                        values_to = "corr") %>%
    dplyr::mutate(sim_var = factor(sim_var, levels = order)) %>%
    dplyr::mutate(obs_var = factor(obs_var, levels = order)) %>%
    ggplot2::ggplot(ggplot2::aes(x = obs_var, y = sim_var, fill = corr)) +
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

  pred_plots$heatmap <- corr_htmp
  return(cowplot::plot_grid(plotlist = pred_plots, nrow = 1))
}
