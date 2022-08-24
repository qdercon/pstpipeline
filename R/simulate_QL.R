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
#' @param gain_loss Fit the dual learning rate model?
#' @param test Simulate test choices in addition to training choices?
#' @param affect Simulate subjective affect ratings (uses full passage-of-time
#' model).
#' @param prev_sample An optional previous sample of id numbers (if you wish to
#' simulate data for the same subset of individual parameters across a number
#' of models).
#' @param raw_df Provide the raw data used to fit the data originally, so that
#' subject IDs can be labelled appropriately.
#' @param ... Other arguments which can be used to control the parameters of the
#' gamma/Gaussian distributions from which parameter values are sampled.
#'
#' @returns Simulated training data (and test data relevant) for a random or
#' previously fitted sample of parameter values.
#'
#' @examples
#' train_sim_2a <- simulate_QL(
#'   sample_size = 10,
#'   alpha_pos_dens = c(shape = 2, scale = 0.1), # default
#'   alpha_neg_dens = c(shape = 2, scale = 0.1), # default
#'   beta_dens = c(mean = 3, sd = 1) # default
#' )
#'
#' @export

simulate_QL <- function(summary_df = NULL,
                        sample_size = NULL,
                        gain_loss = TRUE,
                        test = FALSE,
                        affect = FALSE,
                        prev_sample = NULL,
                        raw_df = NULL,
                        ...) {

  # to appease R CMD check
  variable <- . <- subjID <- value <- trial_no <- id_no <- block <-
    hidden_reward <- question_response <- question_type <- NULL

  if ((!is.null(summary_df) | !is.null(raw_df)) & affect)
    stop("Only random samples supported for affect models.")

  if (is.null(summary_df)) {
    if (is.null(sample_size)) sample_size <- 100
    l <- list(...)
    if (is.null(l$alpha_dens)) l$alpha_dens <- c(2, 0.1) # Gamma(alpha, beta)
    if (is.null(l$alpha_pos_dens)) l$alpha_pos_dens <- l$alpha_dens
    if (is.null(l$alpha_neg_dens)) l$alpha_neg_dens <- l$alpha_dens
    if (is.null(l$beta_dens)) l$beta_dens <- c(3, 1) # Normal(mu, sigma)
    ids_sample <- 1:sample_size

    if (gain_loss) {
      pars_df <- tibble::tibble(
        id_no = ids_sample,
        alpha_pos = rgamma(
          sample_size, shape = l$alpha_pos_dens[1], scale = l$alpha_pos_dens[2]
          ),
        alpha_neg = rgamma(
          sample_size, shape = l$alpha_neg_dens[1], scale = l$alpha_neg_dens[2]
          ),
        beta = rnorm(sample_size, mean = l$beta_dens[1], sd = l$beta_dens[2])
      )
    }
    else {
      pars_df <- tibble::tibble(
        id_no = ids_sample,
        alpha = rgamma(
          sample_size, shape = l$alpha_dens[1], scale = l$alpha_dens[2]
          ),
        beta = rnorm(sample_size, mean = l$beta_dens[1], sd = l$beta_dens[2])
      )
    }
    if (affect) {
      if (is.null(l$gamma_dens)) l$gamma_dens <- c(6, 0.07) # gamma(alpha, beta)
      if (is.null(l$w0_dens)) l$w0_dens <- c(0.4, 0.2) # Normal(mu, sigma)
      if (is.null(l$w1_o_dens)) l$w1_o_dens <- c(-0.05, 0.2) # Normal(mu, sigma)
      if (is.null(l$w1_b_dens)) l$w1_b_dens <- c(-0.1, 0.05) # Normal(mu, sigma)
      if (is.null(l$w2_dens)) l$w2_dens <- c(0.02, 0.05, c(-0.05, Inf))
      if (is.null(l$w3_dens)) l$w3_dens <- c(0.02, 0.05, c(-0.05, Inf))
        # trNormal(mu, sigma, range)

      pars_df <-
        dplyr::bind_rows(
          list("happy" = pars_df, "confident" = pars_df, "engaged" = pars_df),
          .id = "adj"
        ) |>
        dplyr::rowwise() |>
        dplyr::mutate(
          gamma = rgamma(1, shape = l$gamma_dens[1], scale = l$gamma_dens[2]),
          gamma = ifelse(gamma >= 1, runif(1, 0, 1), gamma), # unlikely
          w0    = rnorm(1, l$w0_dens[1], l$w0_dens[2]),
          w1_o  = rnorm(1, l$w1_o_dens[1], l$w1_o_dens[2]),
          w1_b  = rnorm(1, l$w1_b_dens[1], l$w1_b_dens[2]),
          w2    = rnormt(1, l$w2_dens[1], l$w2_dens[2], range = l$w2_dens[3]),
          w3    = rnormt(1, l$w3_dens[1], l$w3_dens[2], range = l$w3_dens[3]])
        ) |>
        dplyr::ungroup()
    }
  }
  else {
    pars_df <- summary_df |>
      dplyr::filter(grepl("alpha|beta", variable)) |>
      dplyr::filter(!grepl("_pr", variable)) |>
      dplyr::filter(!grepl("mu_", variable)) |>
      dplyr::select(variable, mean) |>
      dplyr::mutate(
        id_no = as.numeric(sub("\\].*$", "", sub(".*\\[", "", .[["variable"]])))
        ) |>
      dplyr::mutate(variable = sub("\\[.*$", "", .[["variable"]])) |>
      tidyr::pivot_wider(names_from = variable, values_from = mean)

    if (!is.null(prev_sample)) {
      ids_sample <- prev_sample
    }
    else if (!is.null(sample_size)) { # do we want a sample < everyone
      ids_sample <- sample(1:dim(pars_df)[1], sample_size) |> sort()
    }
    else {
      ids_sample <- 1:dim(pars_df)[1]
    }

    if (!is.null(raw_df)) {
      if (test) raw_df <- raw_df[[1]]
      raw_df <- raw_df |>
        dplyr::select(subjID) |>
        dplyr::distinct() |>
        dplyr::mutate(id_no = dplyr::row_number()) |>
        dplyr::inner_join(
          tibble::as_tibble(ids_sample), by = c("id_no" = "value")
          )
    }
  }

  rewards <- function(i) return (
    rbind(data.frame("trial_block" = i, "type" = rep(12, 20),
                     "hidden_reward" = sample(c(rep(0, 4), rep(1, 16)))),
          data.frame("trial_block" = i, "type" = rep(34, 20),
                     "hidden_reward" = sample(c(rep(0, 6), rep(1, 14)))),
          data.frame("trial_block" = i, "type" = rep(56, 20),
                     "hidden_reward" = sample(c(rep(0, 8), rep(1, 12))))
          )
  )
    # function that returns a random series of hidden rewards the same way as in
    # the task (i.e., in each block, there are exactly 16 rewarded "A", 14
    # rewarded "B", and 12 rewarded "C" symbols

  all_res <- data.frame()
  pb = txtProgressBar(min = 0, max = length(ids_sample), initial = 0, style = 3)
  drop_count <- 0

  for (id in seq_along(ids_sample)) {

    setTxtProgressBar(pb, id)

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
      dplyr::filter(id_no == id)

    if (gain_loss) {
      alpha_pos <- indiv_pars$alpha_pos
      alpha_neg <- indiv_pars$alpha_neg
      beta <- indiv_pars$beta
    }
    else {
      alpha <- indiv_pars$alpha
      beta <- indiv_pars$beta
    }

    hidden_rewards <- dplyr::bind_rows(lapply(1:6, FUN = rewards)) |>
        # 6 blocks
      dplyr::group_by(type) |>
      dplyr::mutate(trial_no_group = dplyr::row_number()) |>
      dplyr::ungroup()

    training_results <- data.frame(
      "id_no" = ids_sample[id],
      "type" = ifelse(conds[,1] == "A", 12, ifelse(conds[,1] == "C", 34, 56)),
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
      if (is.null(l$question_order)) {
        adj_order <- c("happy", "confident", "engaged")
      }
      training_results <- training_results |>
        dplyr::mutate(
          trial_time = NA,
          question_type = rep(sample(adj_order), 120),
          question_response = NA
        )

      trial_time <- rep(0, 360)
      block_time <- rep(0, 360)
      ev_vec <- rep(0, 360)
      pe_vec <- rep(0, 360)

      gamma <- indiv_pars$gamma
      w0    <- indiv_pars$w0
      w1_o  <- indiv_pars$w1_o
      w1_b  <- indiv_pars$w1_b
      w2    <- indiv_pars$w2
      w3    <- indiv_pars$w3

    }

    # Initial Q values
    Q <- data.frame("A" = 0, "B" = 0, "C" = 0, "D" = 0, "E" = 0, "F" = 0)

    for (i in 1:360) {
      Q_t <- Q
      names(Q_t) <- c("Q_a", "Q_b", "Q_c", "Q_d", "Q_e", "Q_f")
      # probability to choose "correct" stimulus
      p_t <- (exp(Q[conds[i,1]] * beta)) / (exp(Q[conds[i,1]] * beta) +
                                              (exp(Q[conds[i,2]] * beta)))

      # make choice
      choice <- sample(c(1, 0), 1, prob = c(p_t, 1-p_t))
      choice_idx <- ifelse(choice == 1, 1, 2)

      # rewarded?
      if (affect) {
        reward <- ifelse(choice == training_results$hidden_reward[i], 1, -1)
      }
      else reward <- ifelse(choice == training_results$hidden_reward[i], 1, 0)
        # i.e. if they choose "correctly" and hidden reward == 1 or
        # incorrectly but hidden reward == 0
        # recoded as -1 for affect models to allow for negative EVs

      ev <- Q[conds[i,choice_idx]]
      pe <- reward - Q[conds[i,choice_idx]]

      # update Q values
      if (gain_loss) { # i.e., were they rewarded?
        if (pe >= 0) Q[conds[i,choice_idx]] <- ev + alpha_pos * pe
        else Q[conds[i,choice_idx]] <- ev + alpha_neg * pe
      }
      else Q[conds[i,choice_idx]] <- ev + alpha * pe

      training_results$choice[i] <- choice
      training_results$reward[i] <- reward

      if (affect) {
        ev_vec[i] <- ev
        pe_vec[i] <- pe
        if (i > 1) {
          if ((i-1) %% 60 == 0) {
            trial_time[i] <-
              trial_time[i-1] + (rnormt(1, c(5, Inf), 90, 45) / 3600)
            block_time[i] <- 0
          }
          else {
            elapsed <- rnormt(1, c(3, Inf), 8.5, 4) / 3600
            trial_time[i] <- trial_time[i-1] + elapsed
            block_time[i] <- block_time[i-1] + elapsed
          }
        }

        q <- grep(training_results$question_type[i], adj_order)

        rating <-
          w0[q] +
          w1_o[q] * trial_time[i] +
          w1_b[q] * block_time[i] +
          w2[q] * sum(sapply(1:i, function(j) gamma[q]^(i-j) * ev_vec[[j]])) +
          w3[q] * sum(sapply(1:i, function(j) gamma[q]^(i-j) * pe_vec[[j]]))

        training_results$question_response[i] <- rating * 100
        training_results$trial_time[i]        <- trial_time[i] * 60 ## in mins
      }
    }

    if (affect) {
      if (!(any(training_results$question_response > 100) &
          any(training_results$question_response < 0))) {
        # changing one may change the other, so skip these cases (v. rare)
        if (any(training_results$question_response > 100) |
          any(training_results$question_response < 0)) {
          indiv_aff <- training_results |>
            dplyr::filter(question_response > 100 | question_response < 0) |>
            with(unique(question_type))
          for (adj in indiv_aff) {
            adj_num <- grep(adj, adj_order)
            if (any(training_results$question_response < 0)) {
              diff <- min(training_results$question_response)
            }
            else if (any(training_results$question_response > 100)) {
              diff <- max(training_results$question_response) - 100
            }
            training_results <- training_results |>
              dplyr::rowwise() |>
              dplyr::mutate(
                question_response =
                  ifelse(question_type == adj, question_response - diff,
                         question_response)
              ) |>
              dplyr::ungroup()
            pars_df[pars_df$id_no == id & pars_df$adj == adj,]$w0 <-
              pars_df[pars_df$id_no == id & pars_df$adj == adj,]$w0 - diff/100
          }
        }
      }
      if (any(training_results$question_response > 100) |
          any(training_results$question_response < 0)) {
        training_results <- NULL
        drop_count <- drop_count + 1
      }
    }

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
        p_t_test <- (exp(Q[conds_test[j,1]] * beta)) /
          (exp(Q[conds_test[j,1]] * beta) + (exp(Q[conds_test[j,2]] * beta)))
          # probability to choose "correct" stimulus
        choice_test <- sample(c(1, 0), 1, prob = c(p_t_test, 1-p_t_test))
          # make choice

        type <- as.numeric(
          paste0(
            match(
              tolower(conds_test[j, 1]), letters[1:6]),
            match(
              tolower(conds_test[j, 2]), letters[1:6])
            )
        )

        if (grepl("12|34|56", type)) {
          test_type <- "training"
        }
        else if (type/10 < 2) {
          test_type <- "chooseA"
        }
        else if (type %% 10 == 2) {
          test_type <- "avoidB"
        }
        else {
          test_type <- "novel"
        }

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

  }

  all_res <- tibble::as_tibble(all_res) |>
    dplyr::select(-hidden_reward)

  if (!is.null(raw_df)) {
    all_res <- all_res |>
      dplyr::left_join(raw_df, by = "id_no")
    pars_df <- pars_df |>
      dplyr::inner_join(raw_df, by = "id_no")
  }
  else if (affect) {
    all_res <- all_res |>
      dplyr::mutate(subjID = id_no) |>
      dplyr::rowwise() |>
      dplyr::mutate(trial_no_block = trial_no - (trial_block-1)*60) |>
      dplyr::mutate(
        question =
          ifelse(question_type == adj_order[1], 1,
                 ifelse(question_type == adj_order[2], 2, 3)
          )
      ) |>
      dplyr::mutate(reward = ifelse(reward == 0, -1, reward)) |>
      dplyr::group_by(subjID, trial_block) |>
      dplyr::mutate(block_time = trial_time - min(trial_time)) |>
      dplyr::group_by(subjID, question_type) |>
      dplyr::mutate(trial_no_q = order(trial_no, decreasing = FALSE)) |>
      dplyr::ungroup()

    pars_df <- pars_df |>
      dplyr::inner_join(
        tibble::as_tibble(ids_sample), by = c("id_no" = "value")
      ) |>
      dplyr::filter(id_no %in% unique(all_res$id_no)) |>
        # only relevant if drop_count > 1
      dplyr::mutate(subjID = id_no)
  }
  else {
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

  if (drop_count > 0) {
    message(paste0(drop_count,
                   " sampled dataset(s) dropped due to extreme values."))
  }
  return(ret)
}
