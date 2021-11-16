#' Quantify and plot associations between learning parameters and variables of interest
#'
#' \code{parameter_glm} is a wrapper around [cmdstan_glm()], taking as inputs
#' summary tables and raw data from [fit_learning_model()], and outputting
#' the results of GLMs quantifying the association between the individual-level posterior
#' means of each parameter and the independent variable(s) of interest. Gamma GLMs with
#' log link functions are used for learning rate models, while standard Gaussian models with
#' identity link are used for models with inverse temperature as the response variable.
#'
#' @param summary_df List of [cmdstanr::summary()] outputs for the fit(s) of interest.
#' @param raw_df List of raw data inputs to the above fits (in the same order). Used to
#' correctly link subject IDs to independent variables.
#' @param var_of_interest Variable of interest.
#' @param covariates Vector of covariates to control for in the GLMs.
#' @param interaction Optional variable to interact with the variable of interest. The
#' GLMs will then be run twice with this variable reverse coded the second time to obtain
#' posterior samples for the variable of interest in both groups. This variable must be binary
#' and only 1 interaction is allowed.
#' @param recode_na Some demographic questions were conditional, and so there exist NAs. This
#' argument allows these terms to be recoded as appropriate (in all binary cases, this should
#' be set to 0).
#' @param factor_scores Given the factor scores were derived separately, this argument allows the
#' \code{data.frame} containing the factor scores to be supplied, so that these factors can
#' be included in models.
#' @param ... Other arguments to pass to [cmdstan_glm()] (e.g., to control
#' number of warm-up and sampling iterations). In addition, use \code{cores} to chnage the
#' number of parallel chains to sample from.
#'
#' @return A [posterior::draws_df()].
#'
#' @importFrom magrittr %>%
#' @importFrom rlang !! :=
#' @export

parameter_glm <- function(summary_df = list(),
                          raw_df = list(),
                          var_of_interest,
                          covariates,
                          interaction = NULL,
                          recode_na = NULL,
                          factor_scores = NULL,
                          ...) {

  l <- list(...)
  if (is.null(l$algorithm)) l$algorithm <- "sampling"
  if (is.null(l$iter_warmup)) l$iter_warmup <- 1000
  if (is.null(l$iter_sampling)) l$iter_sampling <- 1000
  if (is.null(l$chains)) l$chains <- 4
  if (is.null(l$refresh)) l$refresh <- 0
  if (is.null(l$cores)) l$cores <- getOption("mc.cores", 4)
  if (is.null(getOption("mc.cores"))) options(mc.cores = l$cores)

  ## to appease R CMD check
  variable <- parameter <- . <- subjID <- value <- trial_no <- id_no <-  NULL

  all_data <- list()
  for (s in seq_along(summary_df)) {
    all_data[[s]] <- make_par_df(raw_df[[s]], summary_df[[s]])
  }
  all_data <- data.table::rbindlist(all_data, use.names = TRUE)
  if (!is.null(factor_scores)) {
    all_data <- all_data %>%
      dplyr::left_join(factor_scores, by = "subjID")
  }

  formula <- paste0("posterior_mean ~ ",
                    paste0(var_of_interest, collapse = "+"), "+",
                    paste0(covariates, collapse = "+"))
  if (!is.null(interaction)) {
    int_term <- rlang::sym(interaction)
    formula <- paste0(formula, "+", var_of_interest, "*", interaction)
    if (!is.null(recode_na)) {
      all_data <- all_data %>%
        dplyr::mutate(!!int_term := ifelse(is.na(!!int_term), recode_na, !!int_term))
    }
    vals <- unique(all_data[[interaction]])
    if (length(vals) > 2) {
      stop("Interaction term is non-binary. Perhaps NAs need to be recoded?")
    }
    all_data_recode <- all_data %>%
      dplyr::mutate(!!int_term := ifelse(!!int_term == vals[1], vals[2], vals[1]))
    par_ls_recode <- list()
  }
  par_ls <- list()

  mod <- lm(rlang::parse_expr(formula), data = all_data[1:10,])
    # used to get the correct names for betas (as ordering may change e.g., with interactions)
  beta_names <- attr(mod$terms , "term.labels")

  for (par in unique(all_data$parameter)) {
    cmdstan_fit <- cmdstan_glm(
      formula = rlang::parse_expr(formula),
      family = family_ch(par),
      data = dplyr::filter(all_data, parameter == par),
      algorithm = l$algorithm, iter_warmup = l$iter_warmup,
      iter_sampling = l$iter_sampling, chains = l$chains,
      refresh = l$refresh
    )
    par_ls[[par]] <- cmdstan_fit$draws(format = "df") %>%
      dplyr::rename_with(.fn = function(n) return(beta_names[as.numeric(gsub("\\D", "", n))]),
                         .cols = tidyselect::starts_with("beta"))
  }
  if (!is.null(interaction)) {
    for (par in unique(all_data_recode$parameter)) {
      cmdstan_fit <- cmdstan_glm(
        formula = rlang::parse_expr(formula),
        family = family_ch(par),
        data = dplyr::filter(all_data_recode, parameter == par),
        algorithm = l$algorithm, iter_warmup = l$iter_warmup,
        iter_sampling = l$iter_sampling, chains = l$chains,
        refresh = l$refresh
      )
      par_ls_recode[[par]] <- cmdstan_fit$draws(format = "df") %>%
        dplyr::rename_with(.fn = function(n) return(beta_names[as.numeric(gsub("\\D", "", n))]),
                           .cols = tidyselect::starts_with("beta")) %>%
        dplyr::mutate(recode = TRUE)
    }
  }

  pars_df <- data.table::rbindlist(par_ls, idcol = "parameter")
  if (!is.null(interaction)) {
    pars_df <- pars_df %>% dplyr::mutate(recode = FALSE)
    pars_df_recode <- data.table::rbindlist(par_ls_recode, idcol = "parameter")
    pars_df <- dplyr::bind_rows(pars_df, pars_df_recode)
  }
  return(pars_df)
}
