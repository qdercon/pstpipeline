#' @noRd
#' @keywords internal

preprocess_func_train <- function(raw_data, general_info) {
  # Currently class(raw_data) == "data.table"

  # Use general_info of raw_data
  subjs   <- general_info$subjs
  n_subj  <- general_info$n_subj
  t_subjs <- general_info$t_subjs
  t_max   <- general_info$t_max

  # Initialize (model-specific) data arrays
  option1 <- array(-1, c(n_subj, t_max))
  option2 <- array(-1, c(n_subj, t_max))
  choice  <- array(-1, c(n_subj, t_max))
  reward  <- array(-1, c(n_subj, t_max))

  # Write from raw_data to the data arrays
  for (i in 1:n_subj) {
    subj <- subjs[i]
    t <- t_subjs[i]
    DT_subj <- raw_data[raw_data$subjID == subj]
    option1[i, 1:t] <- DT_subj$type %/% 10
    option2[i, 1:t] <- DT_subj$type %% 10
    choice[i, 1:t]  <- DT_subj$choice
    reward[i, 1:t]  <- DT_subj$reward
  }

  # Wrap into a list for Stan
  data_list <- list(
    N       = n_subj,
    T       = t_max,
    Tsubj   = t_subjs,
    option1 = option1,
    option2 = option2,
    choice  = choice,
    reward  = reward
  )

  # Returned data_list will directly be passed to Stan
  return(data_list)
}

preprocess_func_test <- function(raw_data_train, raw_data_test, general_info) {
  # Use general_info of raw_data
  subjs     <- general_info$subjs
  n_subj    <- general_info$n_subj
  t_subjs   <- general_info$t_subjs
  t_max     <- general_info$t_max
  t_subjs_t <- general_info$t_subjs_t
  t_max_t   <- general_info$t_max_t

  # Initialize (model-specific) data arrays
  option1   <- array(-1, c(n_subj, t_max))
  option2   <- array(-1, c(n_subj, t_max))
  choice    <- array(-1, c(n_subj, t_max))
  reward    <- array(-1, c(n_subj, t_max))

  option1_t <- array(-1, c(n_subj, t_max_t))
  option2_t <- array(-1, c(n_subj, t_max_t))
  choice_t  <- array(-1, c(n_subj, t_max_t))

  # Write from raw_data to the data arrays
  for (i in 1:n_subj) {
    subj <- subjs[i]
    t <- t_subjs[i]
    t_t <- t_subjs_t[i]

    DT_subj_train     <- raw_data_train[raw_data_train$subjID == subj]
    option1[i, 1:t]   <- DT_subj_train$type %/% 10
    option2[i, 1:t]   <- DT_subj_train$type %% 10
    choice[i, 1:t]    <- DT_subj_train$choice
    reward[i, 1:t]    <- DT_subj_train$reward

    DT_subj_test      <- raw_data_test[raw_data_test$subjID == subj]
    option1_t[i, 1:t_t] <- DT_subj_test$type %/% 10
    option2_t[i, 1:t_t] <- DT_subj_test$type %% 10
    choice_t[i, 1:t_t]  <- DT_subj_test$choice

  }

  # Wrap into a list for Stan
  data_list <- list(
    N         = n_subj,
    T         = t_max,
    T_t       = t_max_t,
    Tsubj     = t_subjs,
    Tsubj_t   = t_subjs_t,
    option1   = option1,
    option2   = option2,
    choice    = choice,
    reward    = reward,
    option1_t = option1_t,
    option2_t = option2_t,
    choice_t  = choice_t
  )

  # Returned data_list will directly be passed to Stan
  return(data_list)
}
    # adapted from hBayesDM