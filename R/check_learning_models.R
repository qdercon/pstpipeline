#' Diagnostic plots and fit metrics for training and test data models
#'
#' This function is called automatically by [pstpipeline::fit_learning_models()] when \code{model_checks = TRUE},
#' but can also be run separately if desired.
#'
#' @param draws Post-warmup draws - either a [posterior::draws_array()], a [posterior::draws_list()], or a vector
#' of file paths to the .csv output files. May also be a [posterior::draws_df()] but chain-by-chain diagnostics
#' will not be possible.
#' @param mean_pars Output a plot of the mean parameters?
#' @param diagnostic_plots Output diagnostic traces and histograms? Requires the \pkg{bayesplot} package.
#'
#' @importFrom magrittr %>%
#' @export

check_learning_models <- function(draws, mean_pars = TRUE, diagnostic_plots = TRUE, pal = NULL, font = "",
                                  font_size = 11) {

  if (grepl("draws", class(draws)[1])) {
    if (!grepl("df", class(draws)[1])) {
      if (grepl("draws_list", class(draws)[1])) draws <- posterior::as_draws_array(draws)
      mu_pars <- draws[,,grepl("mu_alpha|mu_beta", dimnames(draws)$variable)]
      draws_df <- FALSE
    }
    else {
      suppressWarnings(
        mu_pars <- draws %>%
          dplyr::select(tidyselect::starts_with("mu_")) %>%
          dplyr::select(-tidyselect::contains("pr"))
        )
      mu_pars_df <- mu_pars
      draws_df <- TRUE
      warning("Data given as 'draws_df', so chain-by-chain diagnostics will not be possible.")
    }
  }
  else if (grepl(".csv", draws[1])) {
    mu_pars <-
      tryCatch(
        cmdstanr::read_cmdstan_csv(
          draws, variables = c("mu_alpha_pos", "mu_alpha_neg", "mu_beta")
        ),
        error = function(e) {
          return(
            cmdstanr::read_cmdstan_csv(draws, variables = c("mu_alpha", "mu_beta"))
          )
        }
      )[["post_warmup_draws"]]
    draws_df <- FALSE
  }
  else stop("Unrecognised data format, see help file.")

  if (is.null(pal)) pal <- c("#ffc9b5", "#648767", "#b1ddf1", "#95a7ce", "#987284", "#3d5a80")

  ret <- list()

  if (mean_pars) {
    if (!draws_df) {
      mu_pars_df <- suppressWarnings(
        posterior::as_draws_df(mu_pars) %>%
        dplyr::select(tidyselect::starts_with("mu_")) %>%
        dplyr::select(-tidyselect::contains("pr"))
      )
    }
    pars <- names(mu_pars_df)
    dens_plts <- list()
    dens_plot <- function(df, par, nbins, p, col, font, font_size) {
      rnge <- range(df[par])
      bin_wdth <- diff(rnge) / nbins
      alpha <- grepl("alpha", par)
      xlims <-

      plt <- df %>%
        dplyr::select(all_of(par)) %>%
        dplyr::rename(value = 1) %>%
        ggplot2::ggplot(ggplot2::aes(x = value)) +
        ggplot2::geom_histogram(
          ggplot2::aes(y = ..count.., fill = "value"), colour = "black", alpha = 0.65,
                       binwidth = bin_wdth, position = "identity"
          )  +
        ggplot2::geom_line(
          ggplot2::aes(
            y = (..density..*(dim(df)[1]*bin_wdth)), colour = "value"
          ), size = 1, stat = 'density') +
        ggplot2::scale_colour_manual(values = col) +
        ggplot2::scale_fill_manual(values = col) +
        ggplot2::guides(colour = "none", fill = "none") +
        ggplot2::scale_x_continuous(
          name = bquote(
            .(rlang::parse_expr(paste0(strsplit(par, "_")[[1]][2],
                                ifelse(!alpha, "", paste0("[", strsplit(par, "_")[[1]][3], "]"))
                                )
                         )
              )
          )) +
        ggplot2::ylab("Count") +
        cowplot::theme_half_open(
          font_size = font_size,
          font_family = font
        )
      return(plt)
    }

    for (p in seq_along(pars)) {
      dens_plts[[p]] <- dens_plot(mu_pars_df, nbins = 30, pars[p], col = pal[p], font = font,
                                  font_size = font_size)
    }

    ret$mu_par_dens <- cowplot::plot_grid(plotlist = dens_plts, nrow = 1)
  }
  if (diagnostic_plots) {
    if (length(pal) != 6) pal <- c("#ffc9b5", "#648767", "#b1ddf1", "#95a7ce", "#987284", "#3d5a80")
    bayesplot::bayesplot_theme_set(cowplot::theme_half_open())
    bayesplot::color_scheme_set(pal)

    ret$diagnostics <- list()
    ret$diagnostics$trace <- bayesplot::mcmc_trace(mu_pars)
    ret$diagnostics$rank_hist <- bayesplot::mcmc_rank_hist(mu_pars)
  }

  if (length(ret) == 1 & length(ret[[1]]) == 1) return(ret[[1]])
  else return(ret)

}
