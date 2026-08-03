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
#' @param vb_elbo_retry Should the algorithm retry with different random seeds
#' if the initial ELBO is too low? Applies both to the main VB fit (if
#' \code{vb = TRUE}) and to the VB fit used to generate MCMC initial values
#' (if \code{vb = FALSE} and \code{init} is not otherwise supplied). Defaults
#' to \code{vb}.
#' @param vb_elbo_threshold The minimum acceptable initial ELBO when
#' \code{vb_elbo_retry = TRUE}. Defaults to -9999999.
#' @param vb_max_retries Maximum number of seeds to try if
#' \code{vb_elbo_retry = TRUE}. Defaults to 15.
#' @param vb_sd_check If \code{vb = TRUE}, should the algorithm additionally
#' retry with different seeds when the fitted posterior shows implausibly
#' little between-participant variation in one or more parameters (see
#' [check_par_sd()])? Each attempt requires a full VB refit, so this is
#' considerably more expensive than \code{vb_elbo_retry}. Defaults to
#' \code{vb_elbo_retry}.
#' @param vb_sd_threshold Passed to [check_par_sd()] on each attempt. Defaults
#' to \code{0.01}.
#' @param vb_sd_ignore_pars Passed to [check_par_sd()] on each attempt, for
#' parameters known to be poorly identified (e.g., \code{"w1_o"}). Defaults to
#' none.
#' @param vb_sd_max_retries Maximum number of full VB refits to attempt if
#' \code{vb_sd_check = TRUE}, distinct from \code{vb_max_retries} (which
#' governs the cheap initial-ELBO probe only). Defaults to 10.
#' @param init Initial values for MCMC sampling (ignored if \code{vb = TRUE}).
#' One of:
#' \itemize{
#'   \item \code{NULL} (default): build deterministic default initial values
#'   from each parameter's plausible range (see the internal \code{pars}
#'   list) - group-level means are set from the range's mid-point (inverted
#'   through that parameter's link function), group-level SDs to a small
#'   positive value, and every individual-level parameter to the group mean.
#'   \item A single number: use this as the seed for a variational fit whose
#'   posterior means become the initial values (see \code{"vb"} below), for a
#'   reproducible alternative to \code{vb_elbo_retry}'s seed search.
#'   \item \code{"vb"}: generate initial values from a variational fit,
#'   subject to \code{vb_elbo_retry}/\code{vb_elbo_threshold}/
#'   \code{vb_max_retries} as for the main VB fit (\code{vb_sd_check} is not
#'   applied here).
#'   \item Anything else (e.g. a list or function in one of the formats
#'   accepted by [cmdstanr::sample()]): passed through unchanged.
#' }
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
#' Bayesian metric for model comparison; this requires the \pkg{loo} package,
#' and is computed via [calculate_elpd_loo()].
#' @param save_outputs Save the specified outputs to the disk? Will save to
#' \code{out_dir}.
#' @param return_outputs Include the requested \code{outputs} directly in the
#' returned list? Defaults to the opposite of \code{save_outputs}.
#' Set to \code{FALSE} (only together with \code{save_outputs = TRUE}) to
#' avoid holding large fitted objects (draws, raw data, etc.) in memory -
#' e.g., across many sequential fits in the same notebook/session - in which
#' case the returned list contains the file path each output was saved to
#' instead, to be reloaded later with \code{readRDS()} as needed.
#' @param cores Maximum number of chains to run in parallel. Defaults to
#' \code{options(mc.cores = cores)}
#' or 4 if this is not set (this option will then apply for the rest of the
#' session).
#' @param threads_per_chain Number of threads to use for within-chain
#' parallelisation via \code{reduce_sum}. Defaults to 1 (no within-chain
#' threading). Only the models whose Stan code contains a \code{reduce_sum}
#' call (currently all affect models in \code{rl-affect/}) benefit from this;
#' for other models it is silently ignored. When \code{> 1}, the model is
#' compiled with threading support (\code{cpp_options = list(stan_threads =
#' TRUE)}) and the threads are passed to [cmdstanr::sample()] /
#' [cmdstanr::variational()]. Note the total number of threads used is
#' \code{cores * threads_per_chain} for MCMC, so keep this product at or below
#' the number of physical cores.
#' @param grainsize \code{reduce_sum} grainsize, i.e. the recommended number of
#' participants to sum per partial job. Defaults to 1, which lets the scheduler
#' choose automatically. Only relevant when \code{threads_per_chain > 1}.
#' @param prior_only Sample from the prior only (i.e. skip the likelihood), for
#' prior predictive checks? Only supported by models with a \code{prior_only}
#' data flag (currently all affect models in \code{rl-affect/}). Defaults to
#' \code{FALSE}.
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
                               model = c("1a", "2a", "1a2r", "1a1c"),
                               exp_part = c("training", "test"),
                               affect = FALSE,
                               affect_sfx = c(
                                 "3wt", "4wt_trial", "4wt_block", "4wt_time",
                                 "5wt_time", "delta", "delta-signed"
                               ),
                               adj_order = c("happy", "confident", "engaged"),
                               vb = TRUE,
                               vb_elbo_retry = vb,
                               vb_elbo_threshold = -9999999,
                               vb_max_retries = 100,
                               vb_sd_check = vb_elbo_retry,
                               vb_sd_threshold = 0.05,
                               vb_sd_ignore_pars = character(0),
                               vb_sd_max_retries = 10,
                               init = NULL,
                               ppc = vb,
                               par_recovery = FALSE,
                               task_excl = TRUE,
                               accuracy_excl = FALSE,
                               model_checks = !vb,
                               save_model_as = "",
                               out_dir = "outputs/cmdstan",
                               outputs = c("raw_df", "summary", "draws_list"),
                               save_outputs = TRUE,
                               return_outputs = !save_outputs,
                               cores = getOption("mc.cores", 4),
                               threads_per_chain = cores,
                               grainsize = 1,
                               prior_only = FALSE,
                               ...) {

  if (is.null(getOption("mc.cores"))) options(mc.cores = cores)
  model <- match.arg(model)
  exp_part <- match.arg(exp_part)

  if (threads_per_chain < 1 || threads_per_chain %% 1 != 0) {
    stop("threads_per_chain must be a positive integer.")
  }
  if (grainsize < 1 || grainsize %% 1 != 0) {
    stop("grainsize must be a positive integer.")
  }
  if (!save_outputs && !return_outputs) {
    stop(
      strwrap(
        "return_outputs = FALSE requires save_outputs = TRUE, otherwise
        outputs would be neither saved nor returned.", prefix = " ",
        initial = ""
      )
    )
  }
  use_threads <- threads_per_chain > 1 && affect

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
  # used by find_vb_seed() below regardless of vb, for the cheap ELBO probe
  if (is.null(l$algorithm)) l$algorithm <- "meanfield"
  if (vb) {
    if (is.null(l$iter)) l$iter <- 10000
    if (is.null(l$output_samples)) l$output_samples <- 1000
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
    question_type <- reward <- trial_time <- outc_no <- question_response <-
    trials_elapsed <- NULL

  if (affect) aff_mod <- match.arg(affect_sfx)

  if (!par_recovery) {
    if (task_excl || accuracy_excl) {
      ids <- df_all[["ppt_info"]] |>
        dplyr::select(
          subjID, exclusion, final_block_AB, tidyselect::any_of("distanced")
        )
      if (accuracy_excl) ids <- ids |> dplyr::filter(final_block_AB >= 0.6)
      if (task_excl) ids <- ids |> dplyr::filter(exclusion == 0)
      ids <- ids |> dplyr::select(subjID, tidyselect::any_of("distanced"))
    } else {
      ids <- df_all[["training"]] |>
        dplyr::distinct(subjID, tidyselect::any_of("distanced"))
    }

    training_df <- df_all[["training"]] |>
      dplyr::right_join(tibble::as_tibble(ids), by = c("subjID")) |>
      tidyr::drop_na(choice) # remove timed out trials

    if (exp_part == "test") {
      test_df <- df_all[["test"]] |>
        dplyr::right_join(tibble::as_tibble(ids), by = c("subjID")) |>
        tidyr::drop_na(choice) # remove timed out trials

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
            qn_response_prev = dplyr::lag(question_response, default = -1),
            # the most trials that can be elapsed is 5 if trials aren't missed
            trials_elapsed = pmin(outc_no - dplyr::lag(outc_no, default = 0), 5)
          ) |>
          dplyr::ungroup()

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

  ## get info a la hBayesDM
  if (exp_part == "training") {
    DT_trials <- raw_df[, .N, by = "subjID"]
    subjs     <- DT_trials$subjID
    n_subj    <- length(subjs)
    t_subjs   <- DT_trials$N
    t_max     <- max(t_subjs)

    general_info <- list(subjs, n_subj, t_subjs, t_max)
    names(general_info) <- c("subjs", "n_subj", "t_subjs", "t_max")

    if (affect) {
      # get max number of trials_elapsed for each subject, take the minimum
      general_info[["i_max"]] <-
        min(raw_df[, max(trials_elapsed), by = "subjID"][[2]])
    }
  } else if (exp_part == "test") {
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

  data_cmdstan$run_gq <- ifelse(ppc, 1, 0)
  data_cmdstan$prior_only <- ifelse(prior_only, 1, 0)
  data_cmdstan$grainsize <- as.integer(grainsize)

  if (all(outputs == "stan_datalist")) return(data_cmdstan)

  cmdstanr::check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)

  ## write relevant stan model to memory and preprocess data
  pref <- switch(
    model,
    "1a" = "Q",
    "2a" = "gainloss_Q",
    "1a2r" = "Q_dualsens",
    "1a1c" = "Q_collinspers",
    stop("Invalid model specified.")
  )
  suff <- ifelse(
    !affect, exp_part, paste("plus_affect", aff_mod, sep = "_")
  )

  if (affect) {
    model_dir <- ifelse(
      grepl("2a", model),
      "rl-affect/2-alpha",
      "rl-affect/1-alpha"
    )
  } else {
    model_dir <- ifelse(
      grepl("2a", model), "choice_models/2-alpha", "choice_models/1-alpha"
    )
  }

  stan_rel_path <- paste0(
    "extdata/stan_files/", model_dir, "/pst_", pref, "_", suff, ".stan"
  )
  stan_path <- system.file(stan_rel_path, package = "pstpipeline")

  if (identical(stan_path, "")) {
    stop(
      paste0(
        "Stan model file not found at ", stan_rel_path, ". ",
        "Please check whether you have provided the correct affect_sfx."
      )
    )
  }

  stan_model <- cmdstanr::cmdstan_model(
    stan_path,
    cpp_options = if (use_threads) list(stan_threads = TRUE) else list()
  )

  ## fit variational model if relevant
  best_summary <- NULL

  # Cheap seed search: run a very short VB fit to check the initial ELBO,
  # retrying with new random seeds until one clears vb_elbo_threshold. Used
  # both for the main VB fit below (if vb) and for the VB-derived MCMC
  # initial values further down (if !vb and vb_elbo_retry is set) - the
  # between-participant SD retry (vb_sd_check) is only applied to the former.
  find_vb_seed <- function() {
    if (!vb_elbo_retry) {
      return(list(seed = sample.int(.Machine$integer.max, 1), success = TRUE))
    }

    best_seed <- NULL
    best_elbo <- -Inf

    for (attempt in 1:vb_max_retries) {
      test_seed <- sample.int(.Machine$integer.max, 1)

      attempt_out_dir <- file.path(
        tempdir(), paste0("vb_seed_check_", test_seed)
      )
      if (!dir.exists(attempt_out_dir)) {
        dir.create(attempt_out_dir, recursive = TRUE)
      }

      test_fit <- tryCatch(
        {
          stan_model$variational(
            data = data_cmdstan,
            seed = test_seed,
            iter = 100, # Short run just to get initial ELBO
            refresh = 0,
            output_samples = 1,
            save_latent_dynamics = TRUE,
            show_messages = FALSE,
            show_exceptions = FALSE,
            algorithm = "meanfield",
            threads = if (use_threads) threads_per_chain else NULL,
            output_dir = attempt_out_dir
          )
        }, error = function(e) {
          NULL
        }
      )

      if (!is.null(test_fit)) {
        latent_file <- tryCatch(
          test_fit$latent_dynamics_files()[1],
          error = function(e) ""
        )

        parsed_elbo <- extract_elbo_from_file(latent_file, target_iter = 100)

        elbo_for_check <- parsed_elbo$target_elbo

        if (!is.na(elbo_for_check) && elbo_for_check > vb_elbo_threshold) {
          return(list(seed = test_seed, success = TRUE))
        } else if (!is.na(elbo_for_check) && elbo_for_check > best_elbo) {
          best_seed <- test_seed
          best_elbo <- elbo_for_check
        }
      }
    }
    list(seed = best_seed, success = FALSE)
  }

  if (vb) {
    if (vb_elbo_retry || vb_sd_check) {
      best_fit <- NULL
      best_seed <- NULL
      best_n_collapsed <- Inf
      success <- FALSE
      sd_attempts <- if (vb_sd_check) vb_sd_max_retries else 1

      for (sd_attempt in seq_len(sd_attempts)) {
        seed_result <- find_vb_seed()
        if (is.null(seed_result$seed)) next
        if (!seed_result$success) {
          warning(
            strwrap(
              "Could not find a seed with an acceptable initial ELBO within
              vb_max_retries; trying the best seed found...",
              prefix = " ", initial = ""
            )
          )
        }

        candidate_fit <- tryCatch(
          stan_model$variational(
            data = data_cmdstan,
            seed = seed_result$seed,
            iter = l$iter,
            refresh = if (vb_sd_check) 0 else l$refresh,
            output_samples = l$output_samples,
            algorithm = l$algorithm,
            tol_rel_obj = l$tol_rel_obj,
            show_messages = !vb_sd_check,
            show_exceptions = !vb_sd_check,
            threads = if (use_threads) threads_per_chain else NULL,
            output_dir = out_dir
          ),
          error = function(e) NULL
        )
        if (is.null(candidate_fit)) next

        if (vb_sd_check) {
          stan_vars <- candidate_fit$metadata()$stan_variables
          cand_summary <- candidate_fit$summary(
            variables = stan_vars[is_learning_model_var(stan_vars)]
          )
          n_collapsed <- nrow(
            check_par_sd(cand_summary, vb_sd_threshold, vb_sd_ignore_pars)$collapsed
          )
        } else {
          cand_summary <- NULL
          n_collapsed <- 0
        }

        if (n_collapsed == 0) {
          if (!is.null(best_fit)) unlink(best_fit$output_files())
          best_fit <- candidate_fit
          best_seed <- seed_result$seed
          best_summary <- cand_summary
          success <- TRUE
          if (vb_sd_check) {
            message(
              "VB attempt ", sd_attempt, "/", sd_attempts, ": converged, ",
              "all parameter SDs pass (seed = ", seed_result$seed, ")."
            )
          }
          break
        } else if (n_collapsed < best_n_collapsed) {
          if (!is.null(best_fit)) unlink(best_fit$output_files())
          best_n_collapsed <- n_collapsed
          best_fit <- candidate_fit
          best_seed <- seed_result$seed
          best_summary <- cand_summary
        } else {
          unlink(candidate_fit$output_files())
        }
        if (vb_sd_check) {
          message(
            "VB attempt ", sd_attempt, "/", sd_attempts, ": ", n_collapsed,
            " parameter group(s) below vb_sd_threshold (seed = ",
            seed_result$seed, ") - retrying..."
          )
        }
      }

      if (!success) {
        if (is.null(best_fit)) {
          stop("Could not obtain a valid VB fit within the retry budget.")
        }
        warning(
          strwrap(
            paste0(
              "Could not find a VB fit with no near-zero between-participant
              parameter SDs after ", sd_attempts, " attempt(s); using the
              best fit found (", best_n_collapsed, " parameter group(s)
              still flagged; see vb_sd_threshold / vb_sd_ignore_pars)."
            ), prefix = " ", initial = ""
          )
        )
      }
      fit <- best_fit
      l$seed <- best_seed
    } else {
      fit <- stan_model$variational(
        data = data_cmdstan,
        seed = l$seed,
        iter = l$iter,
        refresh = l$refresh,
        output_samples = l$output_samples,
        algorithm = l$algorithm,
        tol_rel_obj = l$tol_rel_obj,
        threads = if (use_threads) threads_per_chain else NULL,
        output_dir = out_dir
      )
    }
  } else {
    gen_init_vb <- function(model, data_list, parameters, affect, seed = NULL) {
      seed_result <- if (!is.null(seed)) {
        list(seed = seed, success = TRUE)
      } else {
        find_vb_seed()
      }
      if (is.null(seed) && vb_elbo_retry && !seed_result$success) {
        warning(
          strwrap(
            "Could not find an initial-value seed with an acceptable initial
            ELBO within vb_max_retries; using the best seed found for initial
            values...", prefix = " ", initial = ""
          )
        )
      }
      fit_vb <- model$variational(
        data = data_list,
        seed = seed_result$seed,
        refresh = l$refresh,
        threads = if (use_threads) threads_per_chain else NULL
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
          ret
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
            mu_gm = as.vector(m_vb[startsWith(names(m_vb), "mu_gm")]),
            sigma_gm = as.vector(m_vb[startsWith(names(m_vb), "sigma_gm")]),
            aff_mu_phi = as.vector(m_vb[startsWith(names(m_vb), "aff_mu_phi")]),
            aff_sigma_phi = as.vector(
              m_vb[startsWith(names(m_vb), "aff_sigma_phi")]
            )
          )

          matrix_pars <- c(
            "w0", "w1_o", "w1_b", "w2", "w3", "w3_pos", "w3_neg",
            "gm", "phi"
          )
          vec_pars <- setdiff(names(parameters), matrix_pars)

          for (p in vec_pars) {
            ret[[paste0(p, "_pr")]] <-
              as.vector(m_vb[startsWith(names(m_vb), paste0(p, "_pr"))])
          }

          for (q in intersect(names(parameters), matrix_pars)) {
            m_vb_tr <- m_vb[
              names(m_vb[startsWith(names(m_vb), paste0(q, "_pr"))])
            ]
            ret[[paste0(q, "_pr")]] <- cbind(
              as.vector(m_vb_tr[endsWith(names(m_vb_tr), ",1]")]),
              as.vector(m_vb_tr[endsWith(names(m_vb_tr), ",2]")]),
              as.vector(m_vb_tr[endsWith(names(m_vb_tr), ",3]")])
            )
          }
          ret
        }
      }
    }
    if (model == "1a") {
      pars <- list(
        "alpha" = c(0, 0.5, 1),
        "beta" = c(0, 1, 10)
      )
    } else if (model == "1a2r") {
      pars <- list(
        "alpha" = c(0, 0.5, 1),
        "rho_pos" = c(0, 1, 10),
        "rho_neg" = c(0, 1, 10)
      )
    } else if (model == "1a1c") {
      pars <- list(
        "alpha" = c(0, 0.5, 1),
        "pers" = c(0, 0.5, 1),
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
      pars[["w0"]] <- c(-1, 0, 1)
      if (!grepl("3wt", affect_sfx)) pars[["w1_o"]] <- c(-4, 0, 4)
      if (grepl("5wt", affect_sfx)) pars[["w1_b"]] <- c(-2, 0, 2)
      if (!grepl("delta", affect_sfx)) {
        pars[["w2"]] <- c(-2, 0, 2)
        pars[["w3"]] <- c(-2, 0, 2)
      } else {
        pars[["w1_o"]] <- c(-4, 0, 4)
        pars[["w2"]] <- c(-2, 0, 2)
        if (grepl("signed", affect_sfx)) {
          pars[["w3_pos"]] <- c(-2, 0, 2)
          pars[["w3_neg"]] <- c(-2, 0, 2)
        } else {
          pars[["w3"]] <- c(-2, 0, 2)
        }
      }
      pars[["gm"]] <- c(0, 0.5, 1)
      pars[["phi"]] <- c(0, 10, 100)
    }

    # Deterministic default init built from the (lower, mid, upper) ranges in
    # pars above: group-level means come from the mid-point (inverted through
    # that parameter's link function), group-level SDs are a small positive
    # value, and every individual is initialised at the group mean.
    gen_init_default <- function(parameters, affect, n_subj) {
      mid <- vapply(parameters, `[`, numeric(1), 2)

      inv_link <- function(p) {
        val <- mid[[p]]
        if (p %in% c("alpha", "alpha_pos", "alpha_neg", "pers", "gm")) {
          stats::qnorm(val)
        } else if (p == "beta") {
          stats::qnorm(val / 10)
        } else if (p %in% c("rho_pos", "rho_neg", "phi")) {
          log(val)
        } else {
          val # w0, w1_o, w1_b, w2, w3, w3_pos, w3_neg: identity link
        }
      }

      ql_pars <- intersect(
        names(parameters),
        c(
          "alpha", "beta", "alpha_pos", "alpha_neg", "rho_pos", "rho_neg",
          "pers"
        )
      )
      wt_pars <- intersect(
        names(parameters),
        c("w0", "w1_o", "w1_b", "w2", "w3", "w3_pos", "w3_neg")
      )
      # For the delta/delta-signed weight hierarchy, w2/w3(/w3_pos/w3_neg) are
      # the *intercept* of a nested non-centered parameterisation, named
      # "<p>_i_pr" rather than "<p>_pr"; the trial-level second stage has no
      # group hyperparameter to invert, so it's left to CmdStan's own default.
      is_delta <- affect && grepl("delta", affect_sfx)
      nested_wt <- c("w2", "w3", "w3_pos", "w3_neg")

      function() {
        ret <- list()

        if (!affect) {
          ret$mu_pr <- vapply(ql_pars, inv_link, numeric(1))
          ret$sigma <- rep(0.2, length(ql_pars))
          for (p in ql_pars) ret[[paste0(p, "_pr")]] <- rep(0, n_subj)
          return(ret)
        }

        ret$mu_ql <- vapply(ql_pars, inv_link, numeric(1))
        ret$sigma_ql <- rep(0.2, length(ql_pars))
        for (p in ql_pars) ret[[paste0(p, "_pr")]] <- rep(0, n_subj)

        if (length(wt_pars) > 0) {
          wt_mid <- vapply(wt_pars, inv_link, numeric(1))
          ret$mu_wt <- matrix(rep(wt_mid, each = 3), nrow = 3)
          ret$sigma_wt <- matrix(0.5, nrow = 3, ncol = length(wt_pars))
          for (p in wt_pars) {
            p_pr <- if (is_delta && p %in% nested_wt) {
              paste0(p, "_i_pr")
            } else {
              paste0(p, "_pr")
            }
            ret[[p_pr]] <- matrix(0, nrow = n_subj, ncol = 3)
          }
        }
        if ("gm" %in% names(parameters)) {
          ret$mu_gm <- rep(inv_link("gm"), 3)
          ret$sigma_gm <- rep(0.5, 3)
          ret$gm_pr <- matrix(0, nrow = n_subj, ncol = 3)
        }
        if ("phi" %in% names(parameters)) {
          ret$aff_mu_phi <- rep(inv_link("phi"), 3)
          ret$aff_sigma_phi <- rep(0.5, 3)
          ret$phi_pr <- matrix(0, nrow = n_subj, ncol = 3)
        }
        ret
      }
    }

    if (is.null(init)) {
      inits <- gen_init_default(parameters = pars, affect = affect, n_subj = n_subj)
    } else if (is.numeric(init) && length(init) == 1) {
      message(
        "Getting initial values from variational inference (seed = ", init,
        ")..."
      )
      inits <- gen_init_vb(
        model = stan_model, data_list = data_cmdstan, parameters = pars,
        affect = affect, seed = init
      )
    } else if (identical(init, "vb")) {
      message("Getting initial values from variational inference...")
      inits <- gen_init_vb(
        model = stan_model, data_list = data_cmdstan, parameters = pars,
        affect = affect
      )
    } else {
      inits <- init
    }
  }

  ## mcmc sample if relevant
  if (!vb) {
    fit <- stan_model$sample(
      data = data_cmdstan,
      seed = l$seed,
      init = inits,
      refresh = l$refresh, # default = 100
      chains = l$chains, # default = 4
      iter_warmup = l$iter_warmup, # default = 1000
      iter_sampling = l$iter_sampling, # default = 1000
      adapt_delta = l$adapt_delta, # default = 0.8
      step_size = l$step_size, # default = 1
      max_treedepth = l$max_treedepth, # default = 10
      threads_per_chain = if (use_threads) threads_per_chain else NULL,
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

  # Saves (if save_outputs) then either returns value itself, or (if
  # return_outputs = FALSE) just the file path it was saved to - to avoid
  # holding large fitted objects in memory across many sequential fits.
  save_or_return <- function(value, suffix) {
    path <- paste0(out_dir, "/", save_model_as, suffix, ".RDS")
    if (save_outputs) saveRDS(value, file = path)
    if (return_outputs) value else path
  }

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
  if (any(outputs == "model_env")) {
    ret$fit <- if (return_outputs) fit else paste0(out_dir, "/", save_model_as, ".RDS")
  }
  if (any(outputs == "summary")) {
    summary_df <- if (!is.null(best_summary)) {
      best_summary
    } else {
      stan_vars <- fit$metadata()$stan_variables
      fit$summary(variables = stan_vars[is_learning_model_var(stan_vars)])
    }
    ret$summary <- save_or_return(summary_df, "_summary")
  }
  if (any(outputs == "draws_list")) {
    # the least memory intensive format to load
    ret$draws_list <- save_or_return(fit$draws(format = "list"), "_draws_list")
  }
  if (any(outputs == "stan_datalist")) {
    ret$stan_datalist <- save_or_return(data_cmdstan, "_stan_datalist")
  }
  if (any(outputs == "raw_df")) {
    ret$raw_df <- save_or_return(raw_df, "_raw_df")
  }
  if (any(outputs == "loo_obj")) {
    ret$loo_obj <- save_or_return(
      calculate_elpd_loo(fit, cores = cores), "_loo_obj"
    )
  }
  if (any(outputs == "diagnostics") && !vb) {
    ret$diagnostics <- save_or_return(
      fit$cmdstan_diagnose(), "_cmdstan_diagnostics"
    )
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
        ".csv"
      )
    )
  }

  if (!return_outputs) {
    rm(fit, data_cmdstan, raw_df)
    if (exists("best_fit", inherits = FALSE)) rm(best_fit)
    gc(verbose = FALSE)
  }
  ret
}
