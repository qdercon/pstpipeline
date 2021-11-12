#' Simulate data from 1-alpha and 2-alpha Q-learning models
#'
#' \code{simulate_QL} is a function to simulate data from 1-alpha and 2-alpha Q-learning models,
#' with an experiment structure identical to that run online. The parameter values can be from (a
#' sample of) those fitted previously to the real data, or can be randomly sampled.
#'
#' @param summary_df [cmdstanr::summary()] output containing posterior mean parameter values for a group
#' of individuals. If left as NULL, a random sample of size \code{sample_size} will be taken from specified
#' distributions.
#' @param sample_size How may sets of parameters to sample; defaults to 100 or the number of individuals in
#' the \code{summary_df}.
#' @param gain_loss Fit the 2-learning-rate or 1-learning-rate model?
#' @param test Simulate test choices in addition to training choices?
#' @param prev_sample An optional previous sample of id numbers (if you wish to simulate data for the same
#' subset of individual parameters across a number of models).
#' @param raw_df Provide the raw data used to fit the data originally, so that subject IDs can be
#' labelled appropriately.
#' @param ... Other arguments which can be used to control the parameters of the gamma/Gaussian distributions
#' from which parameter values are sampled.
#'
#' @return Simulated training data (and test data relevant) for a random or previously fitted sample of
#' parameter values.
#'
#' @importFrom magrittr %>%
#' @export

simulate_QL <- function(summary_df = NULL, sample_size = NULL, gain_loss = TRUE,
                        test = FALSE, prev_sample = NULL, raw_df = NULL, ...) {

  # to appease R CMD check
  variable <- . <- subjID <- value <- trial_no <- id_no <- block <-
    hidden_reward <- NULL

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
        id_no     = ids_sample,
        alpha_pos = rgamma(sample_size, shape = l$alpha_pos_dens[1], scale = l$alpha_pos_dens[2]),
        alpha_neg = rgamma(sample_size, shape = l$alpha_neg_dens[1], scale = l$alpha_neg_dens[2]),
        beta      = rnorm(sample_size, mean = l$beta_dens[1], sd = l$beta_dens[2])
      )
    }
    else {
      pars_df <- tibble::tibble(
        id_no = ids_sample,
        alpha  = rgamma(sample_size, shape = l$alpha_dens[1], scale = l$alpha_dens[2]),
        beta   = rnorm(sample_size, mean = l$beta_dens[1], sd = l$beta_dens[2])
      )
    }
  }
  else {
    pars_df <- summary_df %>%
      dplyr::filter(grepl("alpha|beta", variable)) %>%
      dplyr::filter(!grepl("_pr", variable)) %>%
      dplyr::filter(!grepl("mu_", variable)) %>%
      dplyr::select(variable, mean) %>%
      dplyr::mutate(id_no = as.numeric(sub("\\].*$", "", sub(".*\\[", "", .[["variable"]])))) %>%
      dplyr::mutate(variable = sub("\\[.*$", "", .[["variable"]])) %>%
      tidyr::pivot_wider(names_from = variable, values_from = mean)

    if (!is.null(prev_sample)) {
      ids_sample <- prev_sample
    }
    else if (!is.null(sample_size)) { # do we want a sample < everyone
      ids_sample <- sample(1:dim(pars_df)[1], sample_size) %>% sort()
    }
    else {
      ids_sample <- 1:dim(pars_df)[1]
    }

    if (!is.null(raw_df)) {
      if (test) raw_df <- raw_df[[1]]
      ids_df <- raw_df %>%
        dplyr::select(subjID) %>%
        dplyr::distinct() %>%
        dplyr::mutate(id_no = dplyr::row_number()) %>%
        dplyr::inner_join(tibble::as_tibble(ids_sample), by = c("id_no" = "value"))
    }
  }

  rewards <- function(i) return (
    rbind(data.frame("block" = i, "type" = rep(12, 20), "hidden_reward" = sample(c(rep(0, 4), rep(1, 16)))),
          data.frame("block" = i, "type" = rep(34, 20), "hidden_reward" = sample(c(rep(0, 6), rep(1, 14)))),
          data.frame("block" = i, "type" = rep(56, 20), "hidden_reward" = sample(c(rep(0, 8), rep(1, 12))))
          )
  )
    # function that returns a random series of hidden rewards the same way as in the task
      # i.e. in each block, there are exactly 16 rewarded "A", 14 rewarded "B", and 12 rewarded "C" symbols

  all_res <- data.frame()
  pb = txtProgressBar(min = 0, max = length(ids_sample), initial = 0, style = 3)

  for (id in seq_along(ids_sample)) {

    setTxtProgressBar(pb, id)

    # make random sequence of trials, number of each condition balanced within blocks
    choice1 <- NULL
    for (i in 1:6) {
      choice1 <- c(choice1, sample(rep(c("A", "C", "E"), each = 20), 60))
    }
    choice2 <- ifelse(choice1 == "A", "B", ifelse(choice1 == "C", "D", "F"))
    conds <- data.frame(choice1, choice2)

    if (test) {
      choice1_test <- tibble::as_tibble(sample(rep(c("A", "C", "D", "E", "F"), times = c(20, 16, 4, 12, 8)))) %>%
        dplyr::rename(choice1_test = value) %>%
        dplyr::mutate(trial_no = dplyr::row_number()) %>%
        dplyr::arrange(choice1_test)

      choice2_test <- c(
        sample(rep(c("B", "C", "D", "E", "F"), each = 4)),
        sample(rep(c("B", "D", "E", "F"), each = 4)),
        rep("B", times = 4),
        sample(rep(c("B", "D", "F"), each = 4)),
        sample(rep(c("B", "D"), each = 4))
      )

      conds_test <- cbind(choice1_test, choice2_test) %>%
        dplyr::arrange(trial_no) %>%
        dplyr::select(-trial_no)
    }

    indiv_pars <- pars_df %>%
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

    hidden_rewards <- dplyr::bind_rows(lapply(1:6, FUN = rewards)) %>% # 6 blocks
      dplyr::group_by(type) %>%
      dplyr::mutate(trial_no_group = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::select(-block)

    training_results <- data.frame(
      "id_no" = ids_sample[id],
      "type" = ifelse(conds[,1] == "A", 12, ifelse(conds[,1] == "C", 34, 56)),
      "choice" = rep(NA, 360),
      "reward" = rep(NA, 360)
    )

    training_results <- training_results %>%
      dplyr::group_by(type) %>%
      dplyr::mutate(trial_no_group = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(hidden_rewards, by = c("type", "trial_no_group"))
    rm(hidden_rewards)

    if (test) {
      training_results <- training_results %>%
        dplyr::mutate(exp_part = "training") %>%
        dplyr::mutate(test_type = NA)
    }

    # Initial Q values
    Q <- data.frame("A" = 0, "B" = 0, "C" = 0,
                    "D" = 0, "E" = 0, "F" = 0)

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
      reward <- ifelse(choice == training_results$hidden_reward[i], 1, 0)
        # i.e. if they choose "correctly" and hidden reward == 1 or incorrectly but hidden reward == 0

      # update Q values
      if (gain_loss) {
        if (reward - Q[conds[i,choice_idx]] >= 0) {
            # i.e., were they rewarded?
          Q[conds[i,choice_idx]] <- Q[conds[i,choice_idx]] +
            alpha_pos * (reward - Q[conds[i,choice_idx]])
        } else {
          Q[conds[i,choice_idx]] <- Q[conds[i,choice_idx]] +
            alpha_neg * (reward - Q[conds[i,choice_idx]])
        }
      } else {
        Q[conds[i,choice_idx]] <- Q[conds[i,choice_idx]] +
          alpha * (reward - Q[conds[i,choice_idx]])
          # single learning rate model
      }

      training_results$choice[i] <- choice
      training_results$reward[i] <- reward

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
        p_t_test <- (exp(Q[conds_test[j,1]] * beta)) / (exp(Q[conds_test[j,1]] * beta) +
                                                          (exp(Q[conds_test[j,2]] * beta)))
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

      test_results <- test_results %>%
        dplyr::group_by(type) %>%
        dplyr::mutate(trial_no_group = dplyr::row_number()) %>%
        dplyr::ungroup()
      all_res <- dplyr::bind_rows(all_res, test_results)
    }
  }

  all_res <- tibble::as_tibble(all_res) %>%
    dplyr::select(-hidden_reward)

  if (!is.null(ids_df)) {
    all_res <- all_res %>%
      dplyr::left_join(ids_df, by = "id_no")
    pars_df <- pars_df %>%
      dplyr::inner_join(ids_df, by = "id_no")
  }
  else {
    all_res <- all_res %>%
      dplyr::mutate(subjID = id_no)
    pars_df <- pars_df %>%
      dplyr::inner_join(tibble::as_tibble(ids_sample), by = c("id_no" = "value")) %>%
      dplyr::mutate(subjID = id_no)
  }

  ret <- list()
  ret$sim <- data.table::as.data.table(all_res)
  ret$pars <- data.table::as.data.table(pars_df)
  return(ret)
}
