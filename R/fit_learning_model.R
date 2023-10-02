#' General function to run Bayesian models using cmdstanr
#'
#' \code{fit_learning_model} uses the package \pkg{cmdstanr}, which is a
#' lightweight R interface to CmdStan. Please note that while it checks if the
#' C++ toolchain is correctly configured, running this function will not install
#' CmdStan itself. This may be as simple as running
#' [cmdstanr::install_cmdstan()], but may require some extra effort (e.g.,
#' pointing R to the install location via [cmdstanr::set_cmdstan_path()]) - see
#' the [cmdstanr vignette](https://mc-stan.org/cmdstanr/articles/cmdstanr.html)
#' for more detail.
#'
#' \code{fit_learning_model} heavily leans on various helper functions from the
#' [\pkg{hBayesDM}](https://ccs-lab.github.io/hBayesDM/) package, and is not as
#' flexible; instead it is designed primarily to be less memory-intensive for
#' our specific use-case and provide only relevant output.
#'
#' @param df_all Raw data outputted from [import_multiple()].
#' @param model Learning model to use, choose from \code{1a} or \code{2a}.
#' @param exp_part Fit to \code{training} or \code{test}?
#' @param affect Fit extended Q-learning model with affect ratings?
#' @param affect_sfx String prefix to identify specific affect model, ignored if
#' \code{affect == FALSE}. Defaults to model with trial-wise passage-of-time.
#' @param adj_order Vector of affect adjectives which is used to define their
#' numerical order in the model output.
#' @param vb Use variational inference to get the approximate posterior? Default
#' is \code{TRUE} for computational efficiency.
#' @param ppc Generate quantities including mean parameters, log likelihood, and
#' posterior predictions? Intended for use with variational algorithm; for MCMC
#' it is recommended to run the separate [generate_posterior_quantities()]
#' function, as this is far less memory intensive.
#' @param par_recovery Method to fit model to simulated data (i.e., from
#' [simulate_QL()]).
#' @param task_excl Apply task-related exclusion criteria (catch questions,
#' digit span = 0)?
#' @param accuracy_excl Apply accuracy-based exclusion criteria (final block AB
#' accuracy >= 0.6)? This is not recommended and is deprecated.
#' @param model_checks Runs [check_learning_models()], returning plots of the
#' group-level posterior densities for the free parameters, and some visual
#' model checks (traceplots of the chains, and rank histograms). Note the visual
#' checks will only be returned if \code{!vb}, as they are only relevant for
#' MCMC fits, and require the \pkg{bayesplot} package.
#' @param save_model_as Name to give to saved model and used to name the .csv
#' files and outputs. Defaults to the Stan model name.
#' @param out_dir Output directory for model fit environment, plus all specified
#' \code{outputs} if \code{save_outputs = TRUE}.
#' @param outputs Specific outputs to return (and save, if \code{save_outputs}).
#' In addition to the defaults, other options are "model_env" (note this is
#' saved automatically, regardless of \code{save_outputs}), and "loo_obj". The
#' latter includes the theoretical expected log-predictive density (ELPD) for a
#' new dataset, plus the leave-one-out information criterion (LOOIC), a fully
#' Bayesian metric for model comparison; this requires the \pkg{loo} package.
#' @param save_outputs Save the specified outputs to the disk? Will save to
#' \code{out_dir}.
#' @param cores Maximum number of chains to run in parallel. Defaults to
#' \code{options(mc.cores = cores)}
#' or 4 if this is not set (this option will then apply for the rest of the
#' session).
#' @param ... Other arguments passed to [cmdstanr::sample()] and/or
#' [check_learning_models]. See the
#' [CmdStan user guide](https://mc-stan.org/docs/2_28/cmdstan-guide/index.html)
#' for full details and defaults.
#'
#' @returns List containing a [cmdstanr::CmdStanVB] or [cmdstanr::CmdStanMCMC]
#' fit object, plus any other outputs passed to \code{outputs}.
#'
#' @importFrom data.table as.data.table .N
#'
#' @examples \dontrun{
#' # Single learning rate Q-learning model fit to training data with MCMC
#'
#' data(example_data)
#' fit1 <- fit_learning_model(
#'   example_data$nd,
#'   model = "1a",
#'   vb = FALSE,
#'   exp_part = "training",
#'   iter_warmup = 1000, # default
#'   iter_sampling = 1000, # default
#'   chains = 4 # default
#' )
#'
#' # Dual learning rate Q-learning model fit to training plus test data with
#' # variational inference
#'
#' data(example_data)
#' fit2 <- fit_learning_model(
#'   example_data$nd,
#'   model = "2a",
#'   exp_part = "test",
#'   vb = TRUE
#' )
#'
#' # Simplest affect model with three weights, fit with variational inference
#'
#' fit3 <- fit_learning_model(
#'   example_data$nd,
#'   model = "2a",
#'   affect = TRUE,
#'   affect_sfx = "3wt",
#'   exp_part = "training",
#'   algorithm = "fullrank"
#' )
#' }
#'
#' @export

fit_learning_model <- function(df_all,
                               model,
                               exp_part,
                               affect = FALSE,
                               affect_sfx = c("3wt", "4wt_trial", "4wt_block",
                                              "4wt_time", "5wt_time",
                                              "cond_prev"),
                               adj_order = c("happy", "confident", "engaged"),
                               vb = TRUE,
                               ppc = vb,
                               par_recovery = FALSE,
                               task_excl = TRUE,
                               accuracy_excl = FALSE,
                               model_checks = !vb,
                               save_model_as = "",
                               out_dir = "outputs/cmdstan",
                               outputs = c("raw_df", "summary", "draws_list"),
                               save_outputs = TRUE,
                               cores = getOption("mc.cores", 4),
                               ...) {

  if (is.null(getOption("mc.cores"))) options(mc.cores = cores)

  if (exp_part == "test" && affect) {
    stop("Affect models will not work for test data.")
  }
  if (affect && !ppc) {
    warning("Separate posterior predictions after affect models not supported.")
    ppc <- TRUE
  }
  if (ppc && !vb) {
    warning(
      strwrap(
        "Loading posterior predictions following MCMC is memory intensive, and
        may result in crashes", prefix = " ", initial = ""
        )
      )
  }
  if (any(outputs == "diagnostics") && vb) {
    warning("Diagnostics are for MCMC only.")
  }

  out_dir <- file.path(getwd(), out_dir)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  l <- list(...)
  if (vb) {
    if (is.null(l$iter)) l$iter <- 10000
    if (is.null(l$output_samples)) l$output_samples <- 1000
    if (is.null(l$algorithm)) l$algorithm <- "meanfield"
    if (is.null(l$tol_rel_obj)) l$tol_rel_obj <- 0.01
  } else { # clearly nothing is being changed, given here just to show defaults
    if (is.null(l$chains)) l$chains <- 4
      # default (explicitly defined here for file naming)
    if (is.null(l$iter_warmup)) l$iter_warmup <- 1000
    # default (explicitly defined here for file naming)
    if (is.null(l$iter_sampling)) l$iter_sampling <- 1000
      # default (explicitly defined here for file naming)
  }

  if (model_checks) {
    if (is.null(l$font)) l$font <- ""
    if (is.null(l$font_size)) l$font_size <- 11
  }

  ## to appease R CMD check
  subjID <- exclusion <- final_block_AB <- choice <- trial_no <- trial_block <-
    question_type <- reward <- trial_time <- question_response <- outc_no <-
    time_elapsed <- NULL

  if (affect) aff_mod <- match.arg(affect_sfx)

  if (!par_recovery) {
    if (task_excl || accuracy_excl) {
      ids <- df_all[["ppt_info"]] |>
        dplyr::select(
          subjID, exclusion, final_block_AB, tidyselect::any_of("distanced"))
      if (accuracy_excl) ids <- ids |> dplyr::filter(final_block_AB >= 0.6)
      if (task_excl) ids <- ids |> dplyr::filter(exclusion == 0)
      ids <- ids |> dplyr::select(subjID, tidyselect::any_of("distanced"))
    } else {
      ids <- df_all[["training"]] |>
        dplyr::distinct(subjID, tidyselect::any_of("distanced"))
    }

    training_df <- df_all[["training"]] |>
      dplyr::right_join(tibble::as_tibble(ids), by = c("subjID")) |>
      tidyr::drop_na(choice)

    if (exp_part == "test") {
      test_df <- df_all[["test"]] |>
        dplyr::right_join(tibble::as_tibble(ids), by = c("subjID")) |>
        tidyr::drop_na(choice)

      raw_df <- list()
      raw_df$train <- data.table::as.data.table(training_df)
      raw_df$test <- data.table::as.data.table(test_df)
    } else {
      if (!affect) {
        raw_df <- data.table::as.data.table(training_df)
      } else {
        training_df <- training_df |>
          dplyr::mutate(trial_no_block = trial_no - (trial_block - 1) * 60) |>
          dplyr::mutate(
            question = dplyr::case_when(
              question_type == adj_order[1] ~ 1,
              question_type == adj_order[2] ~ 2,
              question_type == adj_order[3] ~ 3,
              .default = -1
            )
          ) |>
          dplyr::mutate(reward = ifelse(reward == 0, -1, reward)) |>
          dplyr::group_by(subjID) |>
          dplyr::mutate(outc_no = order(trial_no, decreasing = FALSE)) |>
          dplyr::group_by(trial_block, .add = TRUE) |>
          dplyr::mutate(block_time = trial_time - min(trial_time)) |>
          dplyr::group_by(subjID, question_type) |>
          dplyr::mutate(
            trial_no_q = order(trial_no, decreasing = FALSE),
            qn_response_prev = dplyr::lag(question_response),
            time_elapsed = trial_time - dplyr::lag(trial_time),
            trials_elapsed = outc_no - dplyr::lag(outc_no)
          ) |>
          dplyr::ungroup()

          if (grepl("cond", aff_mod)) {
            training_df <- training_df |> tidyr::drop_na(time_elapsed)
          } else {
            training_df <- training_df |>
              dplyr::mutate(
                qn_response_prev = -1, time_elapsed = -1, trials_elapsed = -1
              )
          }

          raw_df <- data.table::as.data.table(training_df)
      }
    }
  } else {
    if (exp_part == "training") {
      raw_df <- df_all
    } else {
      raw_df <- list()
      raw_df$train <- df_all |>
        dplyr::filter(exp_part == "training") |>
        dplyr::select(-exp_part)
      raw_df$test <- df_all |>
        dplyr::filter(exp_part == "test") |>
        dplyr::select(-exp_part)
    }
  }

  if (all(outputs == "raw_df")) return(raw_df)

  ## get info

  if (exp_part == "training") {
    DT_trials <- raw_df[, .N, by = "subjID"]
    subjs     <- DT_trials$subjID
    n_subj    <- length(subjs)
    t_subjs   <- DT_trials$N
    t_max     <- max(t_subjs)

    general_info <- list(subjs, n_subj, t_subjs, t_max)
    names(general_info) <- c("subjs", "n_subj", "t_subjs", "t_max")
  } else {
    DT_train  <- raw_df$train[, .N, by = "subjID"]
    DT_test   <- raw_df$test[, .N, by = "subjID"]
    subjs     <- DT_train$subjID
    n_subj    <- length(subjs)
    t_subjs   <- DT_train$N
    t_max     <- max(t_subjs)
    t_subjs_t <- DT_test$N
    t_max_t   <- max(t_subjs_t)

    general_info <- list(subjs, n_subj, t_subjs, t_max, t_subjs_t, t_max_t)
    names(general_info) <-
      c("subjs", "n_subj", "t_subjs", "t_max", "t_subjs_t", "t_max_t")
  }

  if (exp_part == "test") {
    data_cmdstan <-
      preprocess_func_test(raw_df$train, raw_df$test, general_info)
  } else {
    if (affect) data_cmdstan <- preprocess_func_affect(raw_df, general_info)
    else data_cmdstan <- preprocess_func_train(raw_df, general_info)
  }

  if (all(outputs == "stan_datalist")) return(data_cmdstan)

  cmdstanr::check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)

  ## write relevant stan model to memory and preprocess data
  label <- ifelse(
    !affect, exp_part, paste("plus_affect", aff_mod, sep = "_")
  )
  stan_model <-
    cmdstanr::cmdstan_model(
      system.file(
        paste0(
          paste(
            "extdata/stan_files/pst",
            ifelse(model == "2a", "gainloss_Q", "Q"), label, sep = "_"),
          ifelse(ppc, "_ppc.stan", ".stan")
        ),
      package = "pstpipeline"
    )
  )

  ## fit variational model if relevant
  if (vb) {
    fit <- stan_model$variational(
      data = data_cmdstan,
      seed = l$seed,
      iter = l$iter,
      refresh = l$refresh,
      output_samples = l$output_samples,
      algorithm = l$algorithm,
      tol_rel_obj = l$tol_rel_obj,
      output_dir = out_dir
    )
  } else if (is.null(l$init)) {
    message("Getting initial values from variational inference...")
    gen_init_vb <- function(model, data_list, parameters, affect) {
      fit_vb <- model$variational(
        data = data_list,
        refresh = l$refresh
      )
      m_vb <- colMeans(fit_vb$draws(format = "df"))

      if (!affect) {
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
      } else {
        function() {
          ret <- list(
            mu_ql = as.vector(m_vb[startsWith(names(m_vb), "mu_ql")]),
            sigma_ql = as.vector(m_vb[startsWith(names(m_vb), "sigma_ql")]),
            mu_wt = rbind(
              as.vector(m_vb[startsWith(names(m_vb), "mu_wt[1,")]),
              as.vector(m_vb[startsWith(names(m_vb), "mu_wt[2,")]),
              as.vector(m_vb[startsWith(names(m_vb), "mu_wt[3,")])
            ),
            sigma_wt = rbind(
              as.vector(m_vb[startsWith(names(m_vb), "sigma_wt[1,")]),
              as.vector(m_vb[startsWith(names(m_vb), "sigma_wt[2,")]),
              as.vector(m_vb[startsWith(names(m_vb), "sigma_wt[3,")])
            ),
            aff_mu_phi = as.vector(m_vb[startsWith(names(m_vb), "aff_mu_phi")]),
            aff_sigma_phi = as.vector(
              m_vb[startsWith(names(m_vb), "aff_sigma_phi")]
            )
          )

          for (p in names(parameters)[1:3]) {
            ret[[paste0(p, "_pr")]] <-
              as.vector(m_vb[startsWith(names(m_vb), paste0(p, "_pr"))])
          }

          for (q in names(parameters)[-(1:3)]) {
            m_vb_tr <- m_vb[
              names(m_vb[startsWith(names(m_vb), paste0(q, "_pr"))])
            ]
            ret[[paste0(q, "_pr")]] <- cbind(
              as.vector(m_vb_tr[endsWith(names(m_vb_tr), ",1]")]),
              as.vector(m_vb_tr[endsWith(names(m_vb_tr), ",2]")]),
              as.vector(m_vb_tr[endsWith(names(m_vb_tr), ",3]")])
            )
          }
          return(ret)
        }
      }
    }
    if (model == "1a") {
      pars <- list(
        "alpha" = c(0, 0.5, 1),
        "beta" = c(0, 1, 10)
      )
    } else {
      pars <- list(
        "alpha_pos" = c(0, 0.5, 1),
        "alpha_neg" = c(0, 0.5, 1),
        "beta" = c(0, 1, 10)
      )
    }
    if (affect) {
      if (!grepl("cond", affect_sfx)) {
        pars[["w0"]] <- c(-1, 0, 1)
        pars[["gm"]] <- c(0, 0.5, 1)
        if (!grepl("3wt", affect_sfx)) {
          pars[["w1_o"]] <- c(-2, 0, 2)
          if (grepl("5wt", affect_sfx)) pars[["w1_b"]] <- c(-2, 0, 2)
        }
      } else {
        pars[["w1_i"]] <- c(-2, 0, 2)
      }
      pars[["w2"]] <- c(-1, 0, 1)
      pars[["w3"]] <- c(-1, 0, 1)
      pars[["phi"]] <- c(0, 10, 100)
    }

    inits <- gen_init_vb(
      model = stan_model,
      data_list = data_cmdstan,
      parameters = pars,
      affect = affect
    )
  }

  ## mcmc sample if relevant
  if (!vb) {
    fit <- stan_model$sample(
      data = data_cmdstan,
      seed = l$seed,
      init = ifelse(is.null(l$init), inits, l$init),
      refresh = l$refresh, # default = 100
      chains = l$chains, # default = 4
      iter_warmup = l$iter_warmup, # default = 1000
      iter_sampling = l$iter_sampling, # default = 1000
      adapt_delta = l$adapt_delta, # default = 0.8
      step_size = l$step_size, # default = 1
      max_treedepth = l$max_treedepth, # default = 10
      output_dir = out_dir
    )
  }

  if (save_model_as == "") {
    save_model_as <- paste(
      "fit_pst", exp_part, model,
      ifelse(vb, "vb", paste0("mcmc_", l$iter_sampling * l$chains)),
      sep = "_"
    )
  }
  fit$save_object(file = paste0(out_dir, "/", save_model_as, ".RDS"))
  ret <- list()
  if (model_checks) {
    if (vb) {
      ret$mu_par_dens <- check_learning_models(
        fit$draws(format = "list"), diagnostic_plots = FALSE, pal = l$pal,
        font = l$font, font_size = l$font_size
      )
    } else {
      ret$model_checks <- list()
      ret$model_checks <- check_learning_models(
        fit$draws(format = "list"), pal = l$pal, font = l$font,
        font_size = l$font_size
      )
    }
  }
  if (any(outputs == "model_env")) ret$fit <- fit
  if (any(outputs == "summary")) {
    ret$summary <- fit$summary()
    if (save_outputs) {
      saveRDS(ret$summary, file = paste0(
        out_dir, "/", save_model_as, "_summary", ".RDS")
        )
    }
  }
  if (any(outputs == "draws_list")) {
    ret$draws <- fit$draws(format = "list")
      # the least memory intensive format to load
    if (save_outputs) {
      saveRDS(ret$draws, file = paste0(
        out_dir, "/", save_model_as, "_draws_list", ".RDS")
        )
    }
  }
  if (any(outputs == "stan_datalist")) {
    ret$stan_datalist <- data_cmdstan
    if (save_outputs) {
      saveRDS(ret$stan_datalist, file = paste0(
        out_dir, "/", save_model_as, "_stan_datalist", ".RDS")
        )
    }
  }
  if (any(outputs == "raw_df")) {
    ret$raw_df <- raw_df
    if (save_outputs) {
      saveRDS(ret$raw_df, file = paste0(
        out_dir, "/", save_model_as, "_raw_df", ".RDS")
        )
    }
  }
  if (any(outputs == "loo_obj")) {
    if (!vb) {
      ret$loo_obj <- fit$loo(cores = cores, save_psis = TRUE)
    } else {
      ll <- ret$draws[[1]][grep("log_lik", names(ret$draws[[1]]))]
      log_lik_mat <- t(do.call(rbind, ll))

      log_p <- ret$draws[[1]]$lp__
      log_g <- ret$draws[[1]]$lp_approx__

      ret$loo_obj <- loo::loo_approximate_posterior(
        log_lik_mat, log_p, log_g, cores = cores, save_psis = TRUE
      )
    }
    if (save_outputs) {
      saveRDS(ret$loo_obj, file = paste0(
        out_dir, "/", save_model_as, "_loo_obj", ".RDS")
        )
    }
  }
  if (any(outputs == "diagnostics") && !vb) {
    ret$diagnostics <- fit$cmdstan_diagnose()
    if (save_outputs) {
      saveRDS(ret$diagnostics, file = paste0(
        out_dir, "/", save_model_as, "_cmdstan_diagnostics", ".RDS")
        )
    }
  }

  ## rename csv output files for improved clarity
  outnames <- fit$output_files()

  for (output in outnames) {
    chain_no <- strsplit(basename(output), "-")[[1]][3]
    file.rename(
      from = output,
      to = paste0(
        out_dir, "/", save_model_as,
        ifelse(vb, paste0("_", l$output_samples), paste0("_chain_", chain_no)),
        ".csv")
      )
  }
  return(ret)
}
