#' Quantify and plot associations between learning parameters and variables of
#' interest
#'
#' \code{parameter_glm} is a wrapper around [cmdstan_glm()], taking as inputs
#' summary tables and raw data from [fit_learning_model()], and outputting
#' the results of GLMs quantifying the association between the individual-level
#' posterior means of each parameter and the independent variable(s) of
#' interest. Gamma GLMs with log link functions are used for learning rate and
#' decay factor models (i.e., positively skewed and constrained between 0 and 1)
#' while standard Gaussian models with identity link functions are used for
#' models with inverse temperature or weights as the response variable.
#'
#' @param summary_df List of [cmdstanr::summary()] outputs for the fit(s) of
#' interest.
#' @param raw_df List of raw data inputs to the above fits (in the same order).
#' Used to correctly link subject IDs to independent variables.
#' @param var_of_interest Variable of interest.
#' @param covariates Vector of covariates to control for in the GLMs.
#' @param affect_number For affect model fits, specify the number (i.e., 1, 2,
#' or 3) of the affect noun/verb of interest. If affect parameters are found in
#' model summaries, and this is not specified, GLMs will default to standard
#' Q-learning parameters.
#' @param interaction Optional variable to interact with the variable of
#' interest. The GLMs will then be run twice with this variable reverse coded
#' the second time to obtain posterior samples for the variable of interest in
#' both groups. This variable must be binary and only 1 interaction is allowed.
#' @param recode_na Some demographic questions were conditional, and so there
#' exist NAs. This argument allows these terms to be recoded as appropriate
#' (in all binary cases, this should be set to 0).
#' @param extra_data Option to supply a data frame with additional derived
#' quantities (e.g., factor scores). Must include a \code{subjID} column.
#' @param rhat_upper,ess_lower Same as [plot_raincloud()].
#' @param ... Other arguments to pass to [cmdstan_glm()] (e.g., to control
#' number of warm-up and sampling iterations). In addition, use \code{cores} to
#' change the number of parallel chains to sample from.
#'
#' @returns A [posterior::draws_df()].
#'
#' @examples \dontrun{
#' # Comparing parameters across groups
#'
#' data(example_data)
#'
#' fit_nd <- fit_learning_model(
#'   example_data$nd,
#'   model = "2a",
#'   vb = FALSE,
#'   exp_part = "training"
#' )
#' fit_dis <- fit_learning_model(
#'   example_data$dis,
#'   model = "2a",
#'   vb = FALSE,
#'   exp_part = "training"
#' )
#'
#' distanced <- parameter_glm(
#'   summary_df = list(fit_nd$summary, fit_dis$summary),
#'   raw_df = list(fit_nd$raw_df, fit_dis$raw_df),
#'   var_of_interest = "distanced",
#'   covariates = c("age", "sex", "digit_span"),
#'   iter_warmup = 1000, iter_sampling = 1000
#' )
#'
#' # Comparing affect model parameters w.r.t. anxiety/depression factor scores
#' # with interaction on distancing
#'
#' factor_scores <- read.csv("data-raw/gillan_scores.csv")[-1] # from Github
#'
#' fit_affect_nd <- fit_learning_model(
#'   example_data$nd,
#'   model = "2a",
#'   affect = TRUE,
#'   exp_part = "training"
#' )
#' fit_affect_dis <- fit_learning_model(
#'   example_data$dis,
#'   model = "2a",
#'   affect = TRUE,
#'   exp_part = "training"
#' )
#'
#' AD_affect_all <- pstpipeline::parameter_glm(
#'   summary_df = list(fit_affect_nd$summary, fit_affect_dis$summary),
#'   raw_df = list(fit_affect_nd$raw_df, fit_affect_dis$raw_df),
#'   var_of_interest = "AD",
#'   covariates = c("age", "sex", "digit_span"),
#'   interaction = "distanced",
#'   affect_number = 1,
#'   extra_data = factor_scores
#' )
#' }
#'
#' @importFrom rlang !! :=
#' @export

parameter_glm <- function(summary_df = list(),
                          raw_df = list(),
                          var_of_interest,
                          covariates,
                          affect_number = NULL,
                          interaction = NULL,
                          recode_na = NULL,
                          extra_data = NULL,
                          rhat_upper = 1.1,
                          ess_lower = 100,
                          ...) {

  l <- list(...)
  if (is.null(l$algorithm)) l$algorithm <- "sampling"
  if (is.null(l$iter_warmup)) l$iter_warmup <- 1000
  if (is.null(l$iter_sampling)) l$iter_sampling <- 1000
  if (is.null(l$chains)) l$chains <- 4
  if (is.null(l$refresh)) l$refresh <- 0
  if (is.null(l$cores)) l$cores <- getOption("mc.cores", 4)
  if (is.null(getOption("mc.cores"))) options(mc.cores = l$cores)
  if (is.null(l$bsl_trnsfm)) l$bsl_trnsfm <- function(x) x

  ## to appease R CMD check
  parameter <- recode <- aff_num <- NULL

  all_data <- list()
  for (s in seq_along(summary_df)) {
    all_data[[s]] <- make_par_df(
      raw_df[[s]], summary_df[[s]], rhat_upper = rhat_upper,
      ess_lower = ess_lower, bsl_trnsfm = l$bsl_trnsfm
    )
  }
  all_data <- data.table::rbindlist(all_data, use.names = TRUE)
  if (!is.null(extra_data)) {
    all_data <- all_data |>
      dplyr::left_join(extra_data, by = "subjID")
  }
  if (!is.null(recode_na)) {
    all_data <- all_data |>
      dplyr::mutate(
        dplyr::across(
          .cols = c(var_of_interest, interaction, covariates),
          .fns = ~ifelse(is.na(.), recode_na, .)
        )
      )
  }
  if ("aff_num" %in% colnames(all_data)) {
    if (is.null(affect_number)) {
      warning(
        "Affect number not specified, defaulting to Q-learning parameters."
      )
      all_data <- all_data |> dplyr::filter(is.na(aff_num))
    } else {
      all_data <- all_data |> dplyr::filter(aff_num == affect_number)
    }
  }

  formula <- paste0("posterior_mean ~ ",
                    paste0(var_of_interest, collapse = "+"), "+",
                    paste0(covariates, collapse = "+"))
  if (!is.null(interaction)) {
    int_term <- rlang::sym(interaction)
    formula <- paste0(formula, "+", var_of_interest, "*", interaction)
    vals <- unique(all_data[[interaction]])
    if (length(vals) > 2) {
      stop("Interaction term is non-binary. Perhaps NAs need to be recoded?")
    }
    all_data_recode <- all_data |>
      dplyr::mutate(
        !!int_term := ifelse(
          !!int_term == vals[1], vals[2], vals[1]
        )
      )
    par_ls_recode <- list()
  }
  par_ls <- list()

  mod <- lm(rlang::parse_expr(formula), data = all_data)
  # used to get the correct names for betas (as ordering may change e.g.,
  # with interactions)
  if (is.factor(all_data[[var_of_interest]])) {
    cov_names <- setdiff(attr(mod$terms, "term.labels"), var_of_interest)
    beta_names <-
      c(names(mod$coefficients)[grep(var_of_interest, names(mod$coefficients))],
        cov_names)

  } else {
    beta_names <- attr(mod$terms, "term.labels")
  }

  for (par in unique(all_data$parameter)) {
    cmdstan_fit <- cmdstan_glm(
      formula = rlang::parse_expr(formula),
      family = family_ch(par),
      data = dplyr::filter(all_data, parameter == par),
      algorithm = l$algorithm, iter_warmup = l$iter_warmup,
      iter_sampling = l$iter_sampling, chains = l$chains,
      refresh = l$refresh
    )
    par_ls[[par]] <- cmdstan_fit$draws(format = "df") |>
      dplyr::rename_with(
        .fn = function(n) return(beta_names[as.numeric(gsub("\\D", "", n))]),
        .cols = tidyselect::starts_with("beta")
      )
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
      par_ls_recode[[par]] <- cmdstan_fit$draws(format = "df") |>
        dplyr::rename_with(
          .fn = function(n) return(beta_names[as.numeric(gsub("\\D", "", n))]),
          .cols = tidyselect::starts_with("beta")
        ) |>
        dplyr::mutate(recode = TRUE)
    }
  }

  pars_df <- data.table::rbindlist(par_ls, idcol = "parameter")
  if (!is.null(interaction)) {
    int_nm <- rlang::sym(paste0(interaction, "_recode"))
    pars_df <- pars_df |> dplyr::mutate(recode = FALSE)
    pars_df_recode <- data.table::rbindlist(par_ls_recode, idcol = "parameter")
    pars_df <- dplyr::bind_rows(pars_df, pars_df_recode) |>
      dplyr::rename(!!int_nm := recode)
  }
  return(pars_df)
}
