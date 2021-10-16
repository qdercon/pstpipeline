#' General function to run Bayesian models using cmdstanr
#'
#' \code{fit_learning_model} uses the package \code{cmdstanr}, which is a lightweight R
#' interface to CmdStan. Please note that while it checks if the C++ toolchain is correctly
#' configured, running this function will not install CmdStan itself. This may be as simple as
#' running \code{cmdstanr::install_cmdstan()}, but may require some extra effort (e.g., pointing
#' R to the install location via \code{cmdstanr::set_cmd_path()}). This vignette explains
#' installation in more detail: https://mc-stan.org/cmdstanr/articles/cmdstanr.html.
#'
#' \code{fit_learning_model} heavily leans on various helper functions from the \code{hBayesDM}
#' package (https://ccs-lab.github.io/hBayesDM/), and is nowhere near as flexible; instead it is
#' designed primarily to be less memory-intensive for our specific use-case and provide only
#' relevant output.
#'
#' @param df List output from import_multiple.
#' @param model Learning model to use, choose from \code{1a} or \code{2a}.
#' @param vb Use variational inference to get the approximate posterior? Default is \code{TRUE}
#' for computational efficiency.
#' @param stan_dir Directory where the stan files are stored.
#' @param test Fit the model to the test data?
#' @param ppc Generate quantities including mean parameters, log likelihood, and posterior predictions?
#' Intended for use with variational algorithm; for MCMC it is recommended to run the separate generate
#' quantities function as this is far less memory intensive.
#' @param diagnostics Outputs and saves diagnostics from \code{$cmdstan_diagnose}. This will not work on
#' saved model objects, nor will it work if \code{vb = TRUE}.
#' @param task_excl Apply task-related exclusion criteria (catch questions, digit span = 0)?
#' @param accuracy_excl Apply accuracy-based exclusion criteria (final block AB accuracy >= 0.6)?
#' @param save_model_as Name to give to saved model and used to name the .csv files and outputs. Defaults
#' to the Stan model name.
#' @param out_dir Output directory for model fit environment, plus all specified \code{outputs} if
#' \code{save_outputs = TRUE}.
#' @param outputs Specific outputs to save separately (in addition to the model environment itself).
#' @param seed Random seed to pass to CmdStan.
#' @param cores Number of CPU cores to run in parallel (passed to \code{options(mc.cores = cores)}).
#' @param ... Other arguments passed to \code{cmdstanr::model$sample}. See the CmdStan user guide for full
#' details and defaults: \code{https://mc-stan.org/docs/2_28/cmdstan-guide/index.html}.
#'
#' @return
#'
#' @importFrom magrittr %>%
#' @export
#'

fit_learning_model <- function(df_all, model, vb = TRUE, stan_dir = "stan_files/", test = FALSE, ppc = vb,
                               diagnostics = !vb, task_excl = TRUE, accuracy_excl = FALSE, save_model_as = "",
                               out_dir = "cmdstan_output/", outputs = c("summary", "draws_df"), save_outputs = TRUE,
                               seed = 123, cores = 4,...) {

  options(mc.cores = cores)
  .datatable.aware = TRUE

  if (ppc & !vb) warning("Loading posterior predictions from MCMC is memory intensive, and may result in crashes.")
  if (diagnostics & vb) warning("Diagnostics are for MCMC only.")

  l <- list(...)
  if (vb) {
    if (is.null(l$iter)) l$iter <- 10000
    if (is.null(l$output_samples)) l$output_samples <- 1000
  }
  else {
    if (is.null(l$init)) l$init <- NULL
    if (is.null(l$refresh)) l$refresh <- NULL # default = 100
    if (is.null(l$chains)) l$chains <- 4 # default (explicitly defined here for file naming)
    if (is.null(l$iter_warmup)) l$iter_warmup <- NULL # default
    if (is.null(l$iter_sampling)) l$iter_sampling <- 1000 # default (explicitly defined here for file naming)
    if (is.null(l$adapt_delta)) l$adapt_delta <- NULL # default = 0.8
    if (is.null(l$step_size)) l$step_size <- NULL # default = 1
    if (is.null(l$max_treedepth)) l$max_treedepth <- NULL # default = 10
  }

  if (task_excl | accuracy_excl) {
    ids <- df_all[["ppt_info"]] %>%
      dplyr::select(subjID, exclusion, final_block_AB)
    if (accuracy_excl) ids <- ids %>% dplyr::filter(final_block_AB >= 0.6)
    if (task_excl) ids <- ids %>% dplyr::filter(exclusion == 0)
    ids <- ids %>% dplyr::select(subjID)
  } else {
    ids <- unique(df_all[["training"]][["subjID"]])
  }

  training_df <- df_all[["training"]] %>%
    dplyr::right_join(tibble::as_tibble(ids), by = c("subjID")) %>%
    tidyr::drop_na(choice)

  if (test) {
    test_df <- df_all[["test"]] %>%
      dplyr::right_join(tibble::as_tibble(ids), by = c("subjID")) %>%
      tidyr::drop_na(choice)

    raw_df_train <- data.table::as.data.table(training_df)
    raw_df_test <- data.table::as.data.table(test_df)

    rm(training_df, test_df)
  }
  else {
    raw_df <- data.table::as.data.table(training_df)
    rm(training_df)
  }

  ## get info

  if (!test) {
    DT_trials <- raw_df[, .N, by = "subjID"]
    subjs     <- DT_trials$subjID
    n_subj    <- length(subjs)
    t_subjs   <- DT_trials$N
    t_max     <- max(t_subjs)

    general_info <- list(subjs, n_subj, t_subjs, t_max)
    names(general_info) <- c("subjs", "n_subj", "t_subjs", "t_max")
  }
  else {
    DT_train  <- raw_df_train[, .N, by = "subjID"]
    DT_test   <- raw_df_test[, .N, by = "subjID"]
    subjs     <- DT_train$subjID
    n_subj    <- length(subjs)
    t_subjs   <- DT_train$N
    t_max     <- max(t_subjs)
    t_subjs_t <- DT_test$N
    t_max_t   <- max(t_subjs_t)

    general_info <- list(subjs, n_subj, t_subjs, t_max, t_subjs_t, t_max_t)
    names(general_info) <- c("subjs", "n_subj", "t_subjs", "t_max", "t_subjs_t", "t_max_t")
  }

  cmdstanr::check_cmdstan_toolchain(fix = TRUE)

  ## write relevant stan model to memory and preprocess data
  if (!test & !ppc) {
    if (model == "1a") stan_model <- cmdstanr::cmdstan_model(paste0(stan_dir,"pst_Q.stan"))
    else if (model == "2a") stan_model <- cmdstanr::cmdstan_model(paste0(stan_dir, "pst_gainloss_Q.stan"))
    data_cmdstan <- pstpipeline::preprocess_func_train(raw_df, general_info)
  } else if (test & !ppc) {
    if (model == "1a") stan_model <- cmdstanr::cmdstan_model(paste0(stan_dir,"pst_Q_test.stan"))
    else if (model == "2a") stan_model <- cmdstanr::cmdstan_model(paste0(stan_dir, "pst_gainloss_Q_test.stan"))
    data_cmdstan <- pstpipeline::preprocess_func_test(raw_df_train, raw_df_test, general_info)
  } else if (!test & ppc) {
    if (model == "1a") stan_model <- cmdstanr::cmdstan_model(paste0(stan_dir,"pst_Q_ppc.stan"))
    else if (model == "2a") stan_model <- cmdstanr::cmdstan_model(paste0(stan_dir, "pst_gainloss_Q_ppc.stan"))
    data_cmdstan <- pstpipeline::preprocess_func_train(raw_df, general_info)
  } else {
    if (model == "1a") stan_model <- cmdstanr::cmdstan_model(paste0(stan_dir,"pst_Q_test_ppc.stan"))
    else if (model == "2a") stan_model <- cmdstanr::cmdstan_model(paste0(stan_dir, "pst_gainloss_Q_test_ppc.stan"))
    data_cmdstan <- pstpipeline::preprocess_func_test(raw_df_train, raw_df_test, general_info)
  }

  ## fit variational model if relevant
  if (vb) {
    fit <- stan_model$variational(
      data = data_cmdstan,
      seed = seed,
      iter = l$iter,
      output_samples = l$output_samples,
      output_dir = out_dir
    )
  }
  else if (is.null(l$init)) {
    message("Getting initial values from variational inference...")
    gen_init_vb <- function(model, data_list, parameters) {
      fit_vb <- model$variational(data = data_list)
      m_vb <- colMeans(as_draws_df(fit_vb$draws()))

      function() {
        ret <- list(
          mu_pr = as.vector(m_vb[startsWith(names(m_vb), "mu_pr")]),
          sigma = as.vector(m_vb[startsWith(names(m_vb), "sigma")])
        )

        for (p in names(parameters)) {
          ret[[paste0(p, "_pr")]] <-
            as.vector(m_vb[startsWith(names(m_vb), paste0(p, "_pr"))])
        }

        return(ret)
      }
    }
    if (model == "1a") {
      pars <- list(
        "alpha" = c(0, 0.5, 1),
        "beta" = c(0, 1, 10)
      )
    }
    else {
      pars <- list(
        "alpha_pos" = c(0, 0.5, 1),
        "alpha_neg" = c(0, 0.5, 1),
        "beta" = c(0, 1, 10)
      )
    }
    inits <- gen_init_vb(
      model = stan_model,
      data_list = data_cmdstan,
      parameters = pars
    )
  }

  ## mcmc sample if relevant
  if (!vb) {
    fit <- stan_model$sample(
      data = data_cmdstan,
      seed = l$seed,
      init = ifelse(is.null(l$init), inits, l$init),
      refresh = l$refresh,
      chains = l$chains,
      iter_warmup = l$iter_warmup,
      iter_sampling = l$iter_sampling,
      adapt_delta = l$adapt_delta,
      step_size = l$step_size,
      max_treedepth = l$max_treedepth,
      output_dir = out_dir
    )
  }

  if (save_model_as == "") save_model_as <- fit$metadata()[["model_name"]]
  fit$save_object(
    file = paste0(out_dir, save_model_as,
                  ifelse(vb, "_vb", paste0("_mcmc_", l$iter_sampling * l$chains)), ".RDS")
    )

  ret <- list()
  if (any(outputs == "model_env")) ret$fit <- fit
  if (any(outputs == "summary")) {
    ret$summary <- fit$summary()
    if (save_outputs) {
      saveRDS(ret$summary, file = paste0(out_dir, save_model_as, "_summary", ".RDS"))
    }
  }
  if (any(outputs == "draws_df")) {
    ret$draws <- posterior::as_draws_df(fit$draws())
    if (save_outputs) {
      saveRDS(ret$draws, file = paste0(out_dir, save_model_as, "_draws", ".RDS"))
    }
  }

  ## rename csv output files for improved clarity
  outnames <- fit$output_files()
  model <- ifelse(save_model_as != "", save_model_as)

  for (output in outnames) {
    chain_no <- strsplit(output, "-")[[1]][3]
    file.rename(
      from = output,
      to = paste0(out_dir, save_model_as,
                  ifelse(vb, paste0("_vb_", l$output_samples, ".csv"),
                         paste0("_mcmc_", l$iter_sampling * l$chains, "chain_", chain_no, ".csv")
                         )
                  )
      )
  }

  return(ret)

}
