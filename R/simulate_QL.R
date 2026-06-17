#' Simulate data from single and dual learning rate Q-learning models
#'
#' \code{simulate_QL} is a function to simulate data from 1-alpha and 2-alpha
#' Q-learning models, with an experiment structure identical to that run online.
#' The parameter values can be from (a sample of) those fitted previously to the
#' real data, or can be randomly sampled.
#'
#' @param summary_df [cmdstanr::summary()] output containing posterior mean
#' parameter values for a group of individuals. If left as NULL, a random sample
#' of size \code{sample_size} will be taken from specified distributions.
#' @param sample_size How may sets of parameters to sample; defaults to 100 or
#' the number of individuals in
#' the \code{summary_df}.
#' @param model Learning model to simulate, one of \code{"1a"}, \code{"2a"},
#' \code{"1a2r"}, or \code{"1a1c"}. Defaults to \code{"2a"}.
#' @param test Simulate test choices in addition to training choices?
#' @param affect Simulate subjective affect ratings (uses full passage-of-time
#' model).
#' @param affect_sfx Affect model suffix to simulate. Mirrors
#' [fit_learning_model()] and can be one of \code{"3wt"}, \code{"4wt_trial"},
#' \code{"4wt_block"}, \code{"4wt_time"}, \code{"5wt_time"},
#' \code{"delta"}, or \code{"delta-signed"}. If left as \code{NULL}, defaults
#' to \code{"4wt_time"}.
#' @param prev_sample An optional previous sample of id numbers (if you wish to
#' simulate data for the same subset of individual parameters across a number
#' of models).
#' @param raw_df Provide the raw data used to fit the data originally, so that
#' subject IDs can be labelled appropriately.
#' @param ... Other arguments which can be used to control the parameters of the
#' Beta/Gaussian distributions from which parameter values are sampled.
#'
#' @returns Simulated training data (and test data relevant) for a random or
#' previously fitted sample of parameter values.
#'
#' @importFrom stats plogis rbeta
#'
#' @examples
#' train_sim_2a <- simulate_QL(
#'   sample_size = 5,
#'   alpha_pos_dens = c(shape = 2, scale = 0.1), # default
#'   alpha_neg_dens = c(shape = 2, scale = 0.1), # default
#'   beta_dens = c(mean = 3, sd = 1) # default
#' )
#'
#' @export

simulate_QL <- function(summary_df = NULL,
                        sample_size = NULL,
                        model = "2a",
                        test = FALSE,
                        affect = FALSE,
                        affect_sfx = NULL,
                        prev_sample = NULL,
                        raw_df = NULL,
                        ...) {

  # to appease R CMD check
  parameter <- subjID <- value <- trial_no <- id_no <- aff_num <-
    hidden_reward <- question_type <- missing_times <- trial_block <- adj <-
    question_response <- outc_lag <- NULL

  l <- list(...)
  sim_model <- match.arg(model, c("1a", "2a", "1a2r", "1a1c"))

  is_dual_alpha <- sim_model == "2a"
  is_dual_sens <- sim_model == "1a2r"
  is_collins_pers <- sim_model == "1a1c"

  if (affect) {
    affect_mod <- if (is.null(affect_sfx)) {
      "4wt_time"
    } else {
      match.arg(
        affect_sfx,
        c(
          "3wt", "4wt_trial", "4wt_block", "4wt_time", "5wt_time",
          "delta", "delta-signed"
        )
      )
    }

    time_prm <- switch(
      affect_mod,
      "3wt" = "none",
      "4wt_trial" = "overall",
      "4wt_block" = "block",
      "4wt_time" = "overall",
      "5wt_time" = c("overall", "block"),
      "delta" = "overall",
      "delta-signed" = "overall"
    )
    is_delta_model <- affect_mod %in% c("delta", "delta-signed")
    signed_delta_model <- affect_mod == "delta-signed"

    sample_time <- is.null(raw_df)
    if (is.null(l$question_order)) {
      l$question_order <- c("happy", "confident", "engaged")
    }
    if (is.null(l$int_max)) l$int_max <- 5
  }

  if (is.null(summary_df)) {
    if (is.null(sample_size)) sample_size <- 100
    if (is.null(l$alpha_dens)) l$alpha_dens <- c(1.5, 3) # Beta(alpha, beta)
    if (is.null(l$alpha_pos_dens)) l$alpha_pos_dens <- l$alpha_dens
    if (is.null(l$alpha_neg_dens)) l$alpha_neg_dens <- l$alpha_dens
    if (is.null(l$beta_dens)) l$beta_dens <- c(3, 4) # Beta(alpha, beta) * 10
    if (is.null(l$rho_pos_dens)) l$rho_pos_dens <- c(3, 4)
    if (is.null(l$rho_neg_dens)) l$rho_neg_dens <- c(3, 4)
    if (is.null(l$pers_dens)) l$pers_dens <- c(0, 0.5) # Normal(mu, sigma)
    ids_sample <- 1:sample_size

    if (is_dual_alpha) {
      pars_df <- tibble::tibble(
        id_no = ids_sample,
        alpha_pos = rbeta(
          sample_size, l$alpha_pos_dens[1], l$alpha_pos_dens[2]
        ),
        alpha_neg = rbeta(
          sample_size, l$alpha_neg_dens[1], l$alpha_neg_dens[2]
        ),
        beta = rbeta(sample_size, l$beta_dens[1], l$beta_dens[2]) * 10
      )
    } else if (is_dual_sens) {
      pars_df <- tibble::tibble(
        id_no = ids_sample,
        alpha = rbeta(sample_size, l$alpha_dens[1], l$alpha_dens[2]),
        rho_pos = rbeta(sample_size, l$rho_pos_dens[1], l$rho_pos_dens[2]) * 10,
        rho_neg = rbeta(sample_size, l$rho_neg_dens[1], l$rho_neg_dens[2]) * 10
      )
    } else if (is_collins_pers) {
      pars_df <- tibble::tibble(
        id_no = ids_sample,
        alpha = rbeta(sample_size, l$alpha_dens[1], l$alpha_dens[2]),
        pers = rnorm(sample_size, l$pers_dens[1], l$pers_dens[2]),
        beta = rbeta(sample_size, l$beta_dens[1], l$beta_dens[2]) * 10
      )
    } else {
      pars_df <- tibble::tibble(
        id_no = ids_sample,
        alpha = rbeta(sample_size, l$alpha_dens[1], l$alpha_dens[2]),
        beta = rbeta(sample_size, l$beta_dens[1], l$beta_dens[2]) * 10
      )
    }
    if (affect) {
      ids_tb <- tibble::tibble(id_no = ids_sample)

      if (is.null(l$w0_dens)) l$w0_dens <- c(0, 0.5) # Normal(mu, sigma)
      if (is.null(l$w1_o_dens)) l$w1_o_dens <- c(-0.5, 1) # Normal(mu, sigma)
      if (is.null(l$w1_b_dens)) l$w1_b_dens <- c(-0.2, 0.5) # Normal(mu, sigma)
      if (is.null(l$w2_dens)) l$w2_dens <- c(0.2, 0.1)
      if (is.null(l$w3_dens)) l$w3_dens <- c(0.2, 0.1)
      if (is.null(l$w3_pos_dens)) l$w3_pos_dens <- l$w3_dens
      if (is.null(l$w3_neg_dens)) l$w3_neg_dens <- l$w3_dens
      if (is.null(l$gamma_dens)) l$gamma_dens <- c(2, 2) # Beta(alpha, beta)

      wt_ids <- dplyr::bind_rows(
        list("happy" = ids_tb, "confident" = ids_tb, "engaged" = ids_tb),
        .id = "adj"
      )

      if (signed_delta_model) {
        wt_df <- wt_ids |>
          dplyr::rowwise() |>
          dplyr::mutate(
            aff_num = grep(adj, l$question_order),
            w0 = rnorm(1, l$w0_dens[1], l$w0_dens[2]),
            w1_o = rnorm(1, l$w1_o_dens[1], l$w1_o_dens[2]),
            w1_b = rnorm(1, l$w1_b_dens[1], l$w1_b_dens[2]),
            w2 = rnorm(1, l$w2_dens[1], l$w2_dens[2]),
            w3_pos = rnorm(1, l$w3_pos_dens[1], l$w3_pos_dens[2]),
            w3_neg = rnorm(1, l$w3_neg_dens[1], l$w3_neg_dens[2]),
            gamma = rbeta(1, l$gamma_dens[1], l$gamma_dens[2])
          ) |>
          dplyr::ungroup()
      } else {
        wt_df <- wt_ids |>
          dplyr::rowwise() |>
          dplyr::mutate(
            aff_num = grep(adj, l$question_order),
            w0 = rnorm(1, l$w0_dens[1], l$w0_dens[2]),
            w1_o = rnorm(1, l$w1_o_dens[1], l$w1_o_dens[2]),
            w1_b = rnorm(1, l$w1_b_dens[1], l$w1_b_dens[2]),
            w2 = rnorm(1, l$w2_dens[1], l$w2_dens[2]),
            w3 = rnorm(1, l$w3_dens[1], l$w3_dens[2]),
            gamma = rbeta(1, l$gamma_dens[1], l$gamma_dens[2])
          ) |>
          dplyr::ungroup()
      }

      if (is_delta_model) {
        wt_a <- wt_df |> dplyr::select(-w0, -w1_o, -w1_b)
        # duplicate wt_b int_max times, with outc_lag = 0, 1, 2, ..., int_max
        if (signed_delta_model) {
          wt_b <- lapply(1:l$int_max, function(i) wt_a) |>
            dplyr::bind_rows(.id = "outc_lag") |>
            dplyr::rowwise() |>
            dplyr::mutate(
              outc_lag = as.integer(outc_lag),
              # lower expectation w/ higher lag
              w2 = w2 * rbeta(1, 1, outc_lag),
              w3_pos = w3_pos * rbeta(1, 1, outc_lag),
              w3_neg = w3_neg * rbeta(1, 1, outc_lag)
            ) |>
            dplyr::ungroup()
        } else {
          wt_b <- lapply(1:l$int_max, function(i) wt_a) |>
            dplyr::bind_rows(.id = "outc_lag") |>
            dplyr::rowwise() |>
            dplyr::mutate(
              outc_lag = as.integer(outc_lag),
              # lower expectation w/ higher lag
              w2 = w2 * rbeta(1, 1, outc_lag),
              w3 = w3 * rbeta(1, 1, outc_lag)
            ) |>
            dplyr::ungroup()
        }
        wt_df <- dplyr::bind_rows(wt_df, wt_b) |> dplyr::select(-gamma)
      }

      pars_df <- dplyr::bind_rows(pars_df, wt_df)

      if (!("block" %in% time_prm)) {
        pars_df  <- pars_df |> dplyr::select(-w1_b)
      }
      if (!("overall" %in% time_prm)) {
        pars_df <- pars_df |> dplyr::select(-w1_o)
      }
    }
  } else {
    pars_df <- clean_summary(summary_df) |>
      dplyr::mutate(
        adj = ifelse(is.na(aff_num), NA, l$question_order[aff_num])
      ) |>
      tidyr::pivot_wider(names_from = parameter, values_from = mean)

    if (!is.null(prev_sample)) {
      ids_sample <- prev_sample
    } else if (!is.null(sample_size)) { # do we want a sample < everyone
      ids_sample <- sample(1:max(pars_df$id_no), sample_size) |> sort()
    } else {
      ids_sample <- seq_len(dim(pars_df))[1]
    }

    if (!is.null(raw_df)) {
      if (test) raw_df <- raw_df[[1]]
      raw_df <- raw_df |>
        dplyr::mutate(
          id_no = as.integer(factor(subjID, levels = unique(raw_df$subjID)))
        ) |>
        dplyr::inner_join(
          tibble::as_tibble(ids_sample), by = c("id_no" = "value")
        )
    }
  }

  rewards <- function(i) {
    (
      rbind(
        data.frame(
          "trial_block" = i, "type" = rep(12, 20),
          "hidden_reward" = sample(c(rep(0, 4), rep(1, 16)))
        ),
        data.frame(
          "trial_block" = i, "type" = rep(34, 20),
          "hidden_reward" = sample(c(rep(0, 6), rep(1, 14)))
        ),
        data.frame(
          "trial_block" = i, "type" = rep(56, 20),
          "hidden_reward" = sample(c(rep(0, 8), rep(1, 12)))
        )
      )
    )
  }
  # function that returns a random series of hidden rewards the same way as in
  # the task (i.e., in each block, there are exactly 16 rewarded "A", 14
  # rewarded "B", and 12 rewarded "C" symbols

  all_res <- data.frame()
  if (is.null(sample_size)) sample_size <- length(ids_sample)
  pb <- txtProgressBar(min = 0, max = sample_size, initial = 0, style = 3)

  for (id in seq_along(ids_sample)) {
    # make random sequence of trials, each condition balanced within blocks
    choice1 <- NULL
    for (i in 1:6) {
      choice1 <- c(choice1, sample(rep(c("A", "C", "E"), each = 20), 60))
    }
    choice2 <- ifelse(choice1 == "A", "B", ifelse(choice1 == "C", "D", "F"))
    conds <- data.frame(choice1, choice2)

    if (test) {
      choice1_test <- tibble::as_tibble(sample(rep(c("A", "C", "D", "E", "F"),
                                                   times = c(20, 16, 4, 12, 8)))
      ) |>
        dplyr::rename(choice1_test = value) |>
        dplyr::mutate(trial_no = dplyr::row_number()) |>
        dplyr::arrange(choice1_test)

      choice2_test <- c(
        sample(rep(c("B", "C", "D", "E", "F"), each = 4)),
        sample(rep(c("B", "D", "E", "F"), each = 4)),
        rep("B", times = 4),
        sample(rep(c("B", "D", "F"), each = 4)),
        sample(rep(c("B", "D"), each = 4))
      )

      conds_test <- cbind(choice1_test, choice2_test) |>
        dplyr::arrange(trial_no) |>
        dplyr::select(-trial_no)
    }

    indiv_pars <- pars_df |>
      dplyr::filter(id_no == ids_sample[id])

    pers <- 0

    if (is_dual_alpha) {
      alpha_pos <- stats::na.omit(indiv_pars$alpha_pos)
      alpha_neg <- stats::na.omit(indiv_pars$alpha_neg)
      beta      <- stats::na.omit(indiv_pars$beta)
    } else if (is_dual_sens) {
      alpha <- stats::na.omit(indiv_pars$alpha)
      rho_pos <- stats::na.omit(indiv_pars$rho_pos)
      rho_neg <- stats::na.omit(indiv_pars$rho_neg)
      beta <- 1
    } else if (is_collins_pers) {
      alpha <- stats::na.omit(indiv_pars$alpha)
      pers <- stats::na.omit(indiv_pars$pers)
      beta <- stats::na.omit(indiv_pars$beta)
    } else {
      alpha <- stats::na.omit(indiv_pars$alpha)
      beta  <- stats::na.omit(indiv_pars$beta)
    }

    hidden_rewards <- dplyr::bind_rows(lapply(1:6, FUN = rewards)) |>
      # 6 blocks
      dplyr::group_by(type) |>
      dplyr::mutate(trial_no_group = dplyr::row_number()) |>
      dplyr::ungroup()

    training_results <- data.frame(
      "id_no" = ids_sample[id],
      "type" = ifelse(conds[, 1] == "A", 12, ifelse(conds[, 1] == "C", 34, 56)),
      "choice" = rep(NA, 360),
      "reward" = rep(NA, 360)
    )

    training_results <- training_results |>
      dplyr::group_by(type) |>
      dplyr::mutate(trial_no_group = dplyr::row_number()) |>
      dplyr::ungroup() |>
      dplyr::left_join(hidden_rewards, by = c("type", "trial_no_group")) |>
      dplyr::mutate(trial_no = dplyr::row_number())
    rm(hidden_rewards)

    if (test) {
      training_results <- training_results |>
        dplyr::mutate(exp_part = "training", test_type = NA)
    }

    if (affect) {
      training_results <- training_results |>
        dplyr::mutate(
          trial_time = NA,
          question_type = as.vector(
            sapply(1:120, function(i) sample(l$question_order))
          ),
          question_response = NA
        )

      if (is_delta_model) {
        training_results <- training_results |>
          dplyr::group_by(question_type) |>
          dplyr::mutate(
            int_trials = trial_no - dplyr::lag(trial_no, default = 0)
          )
      }

      if (!sample_time) {
        time <- raw_df |>
          dplyr::filter(id_no == ids_sample[id]) |>
          dplyr::select(trial_no, block_time, trial_time)
        missing_times <- grep(FALSE, 1:360 %in% time$trial_no)
      }

      ev_vec <- rep(0, 360)
      pe_vec <- rep(0, 360)
      qn_vec <- rep(0, 360)
      trial_time <- rep(0, 360)
      block_time <- rep(0, 360)

      w0  <- stats::na.omit(indiv_pars$w0)
      if (!is_delta_model) {
        gamma <- stats::na.omit(indiv_pars$gamma)
      }
      if ("w1_o" %in% names(indiv_pars))
        w1_o <- stats::na.omit(indiv_pars$w1_o)
      else
        w1_o <- c(0, 0, 0)
      if ("w1_b" %in% names(indiv_pars))
        w1_b <- stats::na.omit(indiv_pars$w1_b)
      else
        w1_b <- c(0, 0, 0)

      if (!is_delta_model) {
        w2 <- stats::na.omit(indiv_pars$w2)
        w3 <- stats::na.omit(indiv_pars$w3)
        signed_w3 <- all(c("w3_pos", "w3_neg") %in% names(indiv_pars))
        if (signed_w3) {
          w3_pos <- stats::na.omit(indiv_pars$w3_pos)
          w3_neg <- stats::na.omit(indiv_pars$w3_neg)
        }
      } else {
        w2 <- sapply(
          1:l$int_max,
          function(i) stats::na.omit(indiv_pars[indiv_pars$outc_lag == i, ]$w2)
        )
        signed_w3 <- all(c("w3_pos", "w3_neg") %in% names(indiv_pars))
        if (signed_w3) {
          w3_pos <- sapply(
            1:l$int_max,
            function(i) {
              stats::na.omit(indiv_pars[indiv_pars$outc_lag == i, ]$w3_pos)
            }
          )
          w3_neg <- sapply(
            1:l$int_max,
            function(i) {
              stats::na.omit(indiv_pars[indiv_pars$outc_lag == i, ]$w3_neg)
            }
          )
        } else {
          w3 <- sapply(
            1:l$int_max,
            function(i) {
              stats::na.omit(indiv_pars[indiv_pars$outc_lag == i, ]$w3)
            }
          )
        }

        int_trials <- training_results$int_trials
      }
    }

    # Initial Q values
    Q <- data.frame("A" = 0, "B" = 0, "C" = 0, "D" = 0, "E" = 0, "F" = 0)
    prev_choice_symbol <- NA_character_

    for (i in 1:360) {

      if (affect) {
        if (!sample_time && i %in% missing_times) next
      }

      val_1 <- Q[[conds[i, 1]]]
      val_2 <- Q[[conds[i, 2]]]
      if (is_dual_sens) {
        val_1 <- ifelse(val_1 >= 0, rho_pos * val_1, rho_neg * val_1)
        val_2 <- ifelse(val_2 >= 0, rho_pos * val_2, rho_neg * val_2)
      }
      rep_1 <- ifelse(is_collins_pers && !is.na(prev_choice_symbol) &&
                        conds[i, 1] == prev_choice_symbol, 1, 0)
      rep_2 <- ifelse(is_collins_pers && !is.na(prev_choice_symbol) &&
                        conds[i, 2] == prev_choice_symbol, 1, 0)

      # probability to choose option1
      p_t <-
        (exp(beta * val_1 + pers * rep_1)) /
        (exp(beta * val_1 + pers * rep_1) + exp(beta * val_2 + pers * rep_2))

      # make choice
      choice <- sample(c(1, 0), 1, prob = c(p_t, 1 - p_t))
      choice_idx <- ifelse(choice == 1, 1, 2)
      chosen_symbol <- as.character(conds[i, choice_idx])

      # rewarded?
      if (affect) {
        reward <- ifelse(choice == training_results$hidden_reward[i], 1, -1)
      } else {
        reward <- ifelse(choice == training_results$hidden_reward[i], 1, 0)
      }
      # i.e. if they choose "correctly" and hidden reward == 1 or
      # incorrectly but hidden reward == 0
      # recoded as -1 for affect models to allow for negative EVs

      ev <- Q[[chosen_symbol]]
      if (is_dual_sens) {
        reward_tr <- ifelse(reward >= 0, rho_pos * reward, rho_neg * reward)
        pe <- reward_tr - Q[[chosen_symbol]]
      } else {
        pe <- reward - Q[[chosen_symbol]]
      }

      # update Q values
      if (is_dual_alpha) {
        if (pe >= 0) {
          Q[[chosen_symbol]] <- ev + alpha_pos * pe
        } else {
          Q[[chosen_symbol]] <- ev + alpha_neg * pe
        }
      } else {
        Q[[chosen_symbol]] <- ev + alpha * pe
      }

      prev_choice_symbol <- chosen_symbol

      training_results$choice[i] <- choice
      training_results$reward[i] <- reward

      if (affect) {
        ev_vec[i] <- ev
        pe_vec[i] <- pe

        if (i > 1 && sample_time) {
          if ((i - 1) %% 60 != 0) {
            elapsed <- 5 + rgamma(1, shape = 3, scale = 2) # after trial
            trial_time[i] <- trial_time[i - 1] + (elapsed / 3600)
            block_time[i] <- 0
          } else {
            elapsed <- 5 + rgamma(1, shape = 5, scale = 5) # after block
            trial_time[i] <- trial_time[i - 1] + (elapsed / 3600)
            block_time[i] <- block_time[i - 1] + (elapsed / 3600)
          }
        } else if (!sample_time) {
          trial_time[i] <- time[time$trial_no == i, ]$trial_time / 60
          block_time[i] <- time[time$trial_no == i, ]$block_time / 60
        }

        q <- grep(training_results$question_type[i], l$question_order)
        qn_vec[i] <- q

        rating <-
          w0[q] +
          w1_o[q] * trial_time[i] +
          w1_b[q] * block_time[i]

        if (!is_delta_model) {
          ev_sum <- sum(
            sapply(1:i, function(j) gamma[q]^(i - j) * ev_vec[[j]])
          )
          if (signed_w3) {
            pe_term <- sum(
              sapply(
                1:i,
                function(j) {
                  w3_j <- ifelse(pe_vec[[j]] >= 0, w3_pos[q], w3_neg[q])
                  gamma[q]^(i - j) * w3_j * pe_vec[[j]]
                }
              )
            )
          } else {
            pe_term <- w3[q] * sum(
              sapply(1:i, function(j) gamma[q]^(i - j) * pe_vec[[j]])
            )
          }

          rating <- rating +
            w2[q] * ev_sum +
            pe_term
        } else {
          for (j in 1:int_trials[i]) {
            t1 <- i + 1 # to include the current trial
            if (signed_w3) {
              w3_t <- ifelse(
                pe_vec[[t1 - j]] >= 0,
                w3_pos[[qn_vec[[t1 - j]], j]],
                w3_neg[[qn_vec[[t1 - j]], j]]
              )
            } else {
              w3_t <- w3[[qn_vec[[t1 - j]], j]]
            }
            rating <- rating +
              w2[[qn_vec[[t1 - j]], j]] * ev_vec[[t1 - j]] +
              w3_t * pe_vec[[t1 - j]]
          }
        }
        training_results$question_response[i] <- plogis(rating) * 100
        training_results$trial_time[i]        <- trial_time[i] * 60 ## in mins
      }
    }

    training_results <- training_results |> tidyr::drop_na(choice)
    all_res <- dplyr::bind_rows(all_res, training_results)

    if (test) {
      test_results <- data.frame(
        "id_no" = ids_sample[id],
        "type" = rep(NA, 60),
        "choice" = rep(NA, 60),
        "reward" = rep(NA, 60),
        "hidden_reward" = rep(NA, 60),
        "exp_part" = rep("test", 60),
        "test_type" = rep(NA, 60)
      )

      for (j in 1:60) {
        val_1_test <- Q[[conds_test[j, 1]]]
        val_2_test <- Q[[conds_test[j, 2]]]
        if (is_dual_sens) {
          val_1_test <- ifelse(
            val_1_test >= 0, rho_pos * val_1_test, rho_neg * val_1_test
          )
          val_2_test <- ifelse(
            val_2_test >= 0, rho_pos * val_2_test, rho_neg * val_2_test
          )
        }
        rep_1_test <- ifelse(
          is_collins_pers && !is.na(prev_choice_symbol) &&
            conds_test[j, 1] == prev_choice_symbol,
          1,
          0
        )
        rep_2_test <- ifelse(
          is_collins_pers && !is.na(prev_choice_symbol) &&
            conds_test[j, 2] == prev_choice_symbol,
          1,
          0
        )
        p_t_test <- (exp(beta * val_1_test + pers * rep_1_test)) /
          (exp(beta * val_1_test + pers * rep_1_test) +
             exp(beta * val_2_test + pers * rep_2_test))
        # probability to choose "correct" stimulus
        choice_test <- sample(c(1, 0), 1, prob = c(p_t_test, 1 - p_t_test))
        # make choice

        chosen_symbol_test <- as.character(
          conds_test[j, ifelse(choice_test == 1, 1, 2)]
        )
        prev_choice_symbol <- chosen_symbol_test

        type <- as.numeric(
          paste0(
            match(tolower(conds_test[j, 1]), letters[1:6]),
            match(tolower(conds_test[j, 2]), letters[1:6])
          )
        )

        if (grepl("12|34|56", type)) test_type <- "training"
        else if (type / 10 < 2) test_type <- "chooseA"
        else if (type %% 10 == 2) test_type <- "avoidB"
        else test_type <- "novel"

        test_results$type[j] <- type
        test_results$choice[j] <- choice_test
        test_results$test_type[j] <- test_type
      }

      test_results <- test_results |>
        dplyr::group_by(type) |>
        dplyr::mutate(trial_no_group = dplyr::row_number()) |>
        dplyr::ungroup()
      all_res <- dplyr::bind_rows(all_res, test_results)
    }
    setTxtProgressBar(pb, id)
  }

  all_res <- tibble::as_tibble(all_res) |>
    dplyr::select(-hidden_reward)

  if (affect) {
    if (!sample_time) {
      pars_df <- pars_df |>
        dplyr::inner_join(
          raw_df |> dplyr::distinct(subjID, id_no), by = "id_no"
        ) |>
        dplyr::mutate(
          id_no = as.integer(factor(id_no, levels = unique(all_res$id_no)))
        )
      # to match up with summary for plotting

      all_res <- all_res |>
        dplyr::left_join(raw_df |> dplyr::distinct(subjID, id_no), by = "id_no")
    } else {
      all_res <- all_res |> dplyr::mutate(subjID = id_no)
    }

    all_res <- all_res |>
      dplyr::rowwise() |>
      dplyr::mutate(trial_no_block = trial_no - (trial_block - 1) * 60) |>
      dplyr::mutate(
        question = ifelse(
          question_type == l$question_order[1], 1,
          ifelse(question_type == l$question_order[2], 2, 3)
        )
      ) |>
      dplyr::mutate(reward = ifelse(reward == 0, -1, reward)) |>
      dplyr::group_by(subjID, trial_block) |>
      dplyr::mutate(block_time = trial_time - min(trial_time)) |>
      dplyr::group_by(subjID, question_type) |>
      dplyr::mutate(
        trial_no_q = order(trial_no, decreasing = FALSE),
        qn_response_prev = dplyr::lag(question_response),
        trials_elapsed = trial_no - dplyr::lag(trial_no, default = 0)
      ) |>
      dplyr::ungroup() |>
      dplyr::arrange(id_no)
  } else if (!is.null(raw_df)) {
    all_res <- all_res |>
      dplyr::left_join(raw_df |> dplyr::distinct(subjID, id_no), by = "id_no")
    pars_df <- pars_df |>
      dplyr::inner_join(raw_df |> dplyr::distinct(subjID, id_no), by = "id_no")
  } else {
    all_res <- all_res |>
      dplyr::mutate(subjID = id_no)
    pars_df <- pars_df |>
      dplyr::inner_join(
        tibble::as_tibble(ids_sample), by = c("id_no" = "value")
      ) |>
      dplyr::mutate(subjID = id_no)
  }

  ret <- list()
  ret$sim <- data.table::as.data.table(all_res)
  ret$pars <- data.table::as.data.table(pars_df)

  ret
}
