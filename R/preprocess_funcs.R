#' @noRd
#' @keywords internal

# adapted from hBayesDM

preprocess_func_train <- function(raw_data,
                                  general_info) {
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

preprocess_func_test <- function(raw_data_train,
                                 raw_data_test,
                                 general_info) {
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

preprocess_func_affect <- function(raw_data, general_info) {
  # Currently class(raw_data) == "data.table"
  # Use general_info of raw_data
  subjs   <- general_info$subjs
  n_subj  <- general_info$n_subj
  t_subjs <- general_info$t_subjs
  t_max   <- general_info$t_max

  # Initialize (model-specific) data arrays
  affect    <- array(0, c(n_subj, t_max))
  afct_prev <- array(0, c(n_subj, t_max))
  block_no  <- array(0, c(n_subj, t_max))
  ovl_trial <- array(0, c(n_subj, t_max))
  int_trial <- array(0, c(n_subj, t_max))
  ovl_time  <- array(0, c(n_subj, t_max))
  blk_time  <- array(0, c(n_subj, t_max))
  int_time  <- array(0, c(n_subj, t_max))
  question  <- array(0, c(n_subj, t_max))

  option1 <- array(-1, c(n_subj, t_max))
  option2 <- array(-1, c(n_subj, t_max))
  choice  <- array(-1, c(n_subj, t_max))
  reward  <- array(-1, c(n_subj, t_max))

  # Write from raw_data to the data arrays
  for (i in 1:n_subj) {
    subj <- subjs[i]
    t <- t_subjs[i]
    DT_subj <- raw_data[raw_data$subjID == subj]
    option1[i, 1:t]   <- DT_subj$type %/% 10
    option2[i, 1:t]   <- DT_subj$type %% 10
    choice[i, 1:t]    <- DT_subj$choice
    reward[i, 1:t]    <- DT_subj$reward
    question[i, 1:t]  <- DT_subj$question
    block_no[i, 1:t]  <- DT_subj$trial_block
    ovl_trial[i, 1:t] <- DT_subj$trial_no_group # ovl question no., by question
    int_trial[i, 1:t] <- DT_subj$trials_elapsed
    ovl_time[i, 1:t]  <- DT_subj$trial_time / 60 # in hours
    blk_time[i, 1:t]  <- DT_subj$block_time / 60
    int_time[i, 1:t]  <- DT_subj$time_elapsed # in mins
    affect[i, 1:t]    <- DT_subj$question_response / 100
    afct_prev[i, 1:t] <- DT_subj$qn_response_prev / 100
  }

  # Wrap into a list for Stan
  data_list <- list(
    N           = n_subj,
    T           = t_max,
    Tsubj       = t_subjs,
    option1     = option1,
    option2     = option2,
    choice      = choice,
    reward      = reward,
    affect      = affect,
    affect_prev = afct_prev, # only relevant for conditional model
    block_no    = block_no,
    ovl_trial   = ovl_trial,
    int_trials  = int_trial, # only relevant for conditional model
    ovl_time    = ovl_time,
    blk_time    = blk_time,
    int_time    = int_time, # only relevant for conditional model
    question    = question
  )

  # Returned data_list will directly be passed to Stan
  return(data_list)
}
