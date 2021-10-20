#' @noRd
#' @keywords internal
#' @importFrom data.table rbindlist
#' @importFrom magrittr %>%
#' @importFrom rlang := !!

import_save_chains <- function(out_files, out_dir, n_indiv, n_trials, prefix, sub_dir = "",
                               vars="y_pred", reorder = TRUE, split_indiv = FALSE) {

  message("Importing chain-by-chain and saving predictions for each individual...")

  if (reorder) {
    var_order <- as.vector(
      sapply(1:n_indiv, function(i)
      sapply(1:n_trials, function(j) paste0("y_pred[", i, ",", j, "]"))
      )
    )
  }

  for (i in seq_along(out_files)) {
    indiv_draws <- cmdstanr::read_cmdstan_csv(out_files[i], variables = vars, format = "draws_list")
    indiv_draws <- indiv_draws$generated_quantities$`1`
    if (reorder) {
      indiv_draws <- indiv_draws %>%
        dplyr::select(tidyselect::all_of(var_order))
    }
    indiv_draws <- split.default(indiv_draws, rep(1:n_indiv, each = n_trials))
    if (split_indiv) {
      lapply(1:n_indiv,
             function(j) saveRDS(indiv_draws[[j]],
                                 paste0(sub_dir, prefix, "_chain_", i, "_ppc_list_", j, ".RDS")))
      rm(indiv_draws)
      gc(reset = TRUE)

      message("\nCompiling individual-level predictions...")

      for (i in 1:n_indiv) {
        dts <- list()
        for (j in seq_along(out_files)) {
          dts[[j]] <- readRDS(paste0(sub_dir, prefix, "_chain_", j, "_ppc_list_", i, ".RDS"))
        }
        indiv_preds <- rbindlist(dts)
        saveRDS(indiv_preds, paste0(sub_dir, prefix, "_ppc_list_all_chains_", i, ".RDS"))
      }
    }
    else {
      saveRDS(indiv_draws, paste0(sub_dir, "/", prefix, "chain_", i, "_ppc_list.RDS"))
    }
  }
}

prep_data <- function(obs, pred_type = c("AB", "CD", "EF"), n_trials, n_draws_chain, prefix, sub_dir,
                      separated_indiv = FALSE, chain = 1, ret_obs = TRUE) {

  indiv <- unique(obs$subjID)
  n_indiv <- length(indiv)

  for (t in seq_along(pred_type)) {
    eval(parse(
      text = paste0(
        paste0("obs_", tolower(pred_type[t])), " <- matrix(NA, nrow = n_indiv, ncol = n_trials/3)"
        )
      ))
    eval(parse(
      text = paste0(
        paste0("pred_", tolower(pred_type[t])), " <- array(dim = c(n_draws_chain, n_trials/3, n_indiv))"
        )
      ))
  }

  if (!separated_indiv) {
    pred_indiv_all <- readRDS(paste0(sub_dir, "/", prefix, "chain_", chain, "_ppc_list.RDS"))
  }

  for (i in 1:n_indiv) {
    if (separated_indiv) {
      pred_indiv <- readRDS(paste0(sub_dir, "/", prefix, "chain_", chain, "_ppc_list_", i, ".RDS"))
      pred_indiv <- pred_indiv[,colSums(pred_indiv) >= 0]
    } else {
      pred_indiv <- pred_indiv_all[[i]]
      pred_indiv <- pred_indiv[,colSums(pred_indiv) >= 0] # columns missing PPCs are set to -1
    }

    missing <- which(!seq(1:n_trials) %in% obs[obs$subjID == indiv[i],]$trial_no)
    for (m in missing) {
      pred_indiv <- pred_indiv %>%
        tibble::add_column(NA, .after = m-1)
    }
    colnames(pred_indiv) <- sapply(1:n_trials, function(n) paste0("y_pred_", n))

    if (any(pred_type == "AB")) {
      ## add observed/predicted choices at correct indexes
      ab_obs <- obs[obs$subjID == indiv[i] & obs$type == 12, ]$choice
      ab_pred_nms <- paste0("y_pred_", obs[obs$subjID == indiv[i] & obs$type == 12,]$trial_no)
      ab_preds <- pred_indiv %>% dplyr::select(tidyselect::all_of(ab_pred_nms))

      ab_trls <- which(seq(1:(n_trials/3)) %in% obs[obs$subjID == indiv[i] & obs$type == 12,]$trial_no_group)
      colnames(ab_preds) <- sapply(ab_trls, function(abt) paste0("y_pred_ab_", abt))

      ab_miss <- which(!seq(1:(n_trials/3)) %in% ab_trls)
      for (abm in ab_miss) {
        ab_obs <- append(ab_obs, as.integer(-1), after = abm-1)
        clmnme <- rlang::sym(paste0("y_pred_ab_", abm))
        ab_preds <- ab_preds %>%
          tibble::add_column(!!clmnme := -1, .after = abm-1)
      }

      obs_ab[i,] <- ab_obs
      pred_ab[,,i] <- as.matrix(ab_preds)
    }

    if (any(pred_type == "CD")) {
      ## add observed/predicted choices at correct indexes
      cd_obs <- obs[obs$subjID == indiv[i] & obs$type == 34, ]$choice
      cd_pred_nms <- paste0("y_pred_", obs[obs$subjID == indiv[i] & obs$type == 34,]$trial_no)
      cd_preds <- pred_indiv %>% dplyr::select(tidyselect::all_of(cd_pred_nms))

      cd_trls <- which(seq(1:(n_trials/3)) %in% obs[obs$subjID == indiv[i] & obs$type == 34,]$trial_no_group)
      colnames(cd_preds) <- sapply(cd_trls, function(cdt) paste0("y_pred_cd_", cdt))

      cd_miss <- which(!seq(1:(n_trials/3)) %in% cd_trls)
      for (cdm in cd_miss) {
        cd_obs <- append(cd_obs, as.integer(-1), after = cdm-1)
        clmnme <- rlang::sym(paste0("y_pred_cd_", cdm))
        cd_preds <- cd_preds %>%
          tibble::add_column(!!clmnme := -1, .after = cdm-1)
      }

      obs_cd[i,] <- cd_obs
      pred_cd[,,i] <- as.matrix(cd_preds)
    }

    if (any(pred_type == "EF")) {
      ef_obs <- obs[obs$subjID == indiv[i] & obs$type == 56, ]$choice
      ef_pred_nms <- paste0("y_pred_", obs[obs$subjID == indiv[i] & obs$type == 56,]$trial_no)
      ef_preds <- pred_indiv %>% dplyr::select(tidyselect::all_of(ef_pred_nms))

      ef_trls <- which(seq(1:(n_trials/3)) %in% obs[obs$subjID == indiv[i] & obs$type == 56,]$trial_no_group)
      colnames(ef_preds) <- sapply(ef_trls, function(eft) paste0("y_pred_ef_", eft))

      ef_miss <- which(!seq(1:(n_trials/3)) %in% ef_trls)
      for (efm in ef_miss) {
        ef_obs <- append(ef_obs, as.integer(-1), after = efm-1)
        clmnme <- rlang::sym(paste0("y_pred_ef_", efm))
        ef_preds <- ef_preds %>%
          tibble::add_column(!!clmnme := -1, .after = efm-1)
      }

      obs_ef[i,] <- ef_obs
      pred_ef[,,i] <- as.matrix(ef_preds)
    }
  }

  obs_list <- list()
  pred_list <- list()

  if (any(pred_type == "AB")) {
    obs_list$obs_AB <- obs_AB
    pred_list$pred_AB <- pred_AB
  }

  if (any(pred_type == "CD")) {
    obs_list$obs_CD <- obs_CD
    pred_list$pred_CD <- pred_CD
  }

  if (any(pred_type == "EF")) {
    obs_list$obs_EF <- obs_EF
    pred_list$pred_EF <- pred_EF
  }

  ret <- list(obs_list, pred_list)
  names(ret) <- c("obs_list", "pred_list")

  return(ret)
}

# prep_avg_preds <- function(preds = NULL, types = c("AB", "CD", "EF"), chains = 4, splits,
#                            n_trials_type = 120, n_subj = 379) {
#
#   if (is.null(preds)) {
#     preds <-
#       sapply(1:chains, function(i) sapply(seq_along(types),
#                                           function (j) paste0("chain", i, "_pred_", types[j])))
#   }
#
#   subj_split = n_trials_type/splits
#     # subjects are averaged across trials
#   trial_split = n_subj/splits
#     # trials are averaged across subjects
#
#   ### SUBJECT-LEVEL
#
#   message("Subject-level averages...")
#
#   pred_subjMCMC <- vector("list", length(types))
#   names(pred_subjMCMC) <- sapply(seq_along(types), function(p) paste0("pred_", types[p], "_subjMCMC"))
#
#   for (t in seq_along(types)) {
#     cat("Trial type ", types[t], "\n")
#
#     eval(parse(text = paste0("pred_subjMCMC_splits <- vector(\"list\", ", splits, ")")))
#
#     for (s in 1:splits) {
#
#       index1_subj = ifelse(s==1, 1, ((s-1)*subj_split)+1)
#       index2_subj = s*subj_split
#
#       subj_list <- vector("list", chains)
#
#       for (c in 1:chains) {
#         eval(parse(text = paste0("pred_chain_", c, " <- readRDS(", "\"", preds[t,c], ".RDS\")")))
#         eval(parse(text = paste0("subj_list[[", c, "]] <- pred_chain_", c, "[,", index1_subj, ":", index2_subj, ",]")))
#         eval(parse(text = paste0("rm(\"pred_chain_", c, "\")")))
#       }
#
#       pred_subj <- abind::abind(subj_list, along=1)
#       rm(subj_list)
#       eval(parse(text = paste0("pred_subjMCMC_splits[[", s, "]] <- apply(pred_subj, c(1,3), mean, na.rm = TRUE)")))
#       rm(pred_subj)
#
#     }
#
#     pred_subjMCMC_splits <- abind::abind(pred_subjMCMC_splits, along=3)
#     pred_subjMCMC[[t]] <- apply(pred_subjMCMC_splits, c(1,2), mean, na.rm = TRUE)
#     rm(pred_subjMCMC_splits)
#
#   }
#
#   ### TRIAL-LEVEL
#
#   message("Trial-level averages...")
#
#   pred_trialMCMC <- vector("list", length(types))
#   names(pred_trialMCMC) <- sapply(seq_along(types), function(p) paste0("pred_", types[p], "_trialMCMC"))
#
#   for (t in seq_along(types)) {
#     cat("Trial type ", types[t], "\n")
#
#     eval(parse(text = paste0("pred_trialMCMC_splits <- vector(\"list\", ", splits, ")")))
#
#     for (s in 1:splits) {
#
#       index1_trial = ifelse(s==1, 1, ((s-1)*trial_split)+1)
#       index2_trial = s*trial_split
#
#       trial_list <- vector("list", chains)
#
#       for (c in 1:chains) {
#         eval(parse(text = paste0("pred_chain_", c, " <- readRDS(", "\"", preds[t,c], ".RDS\")")))
#         eval(parse(text = paste0("trial_list[[", c, "]] <- pred_chain_", c, "[,,", index1_trial, ":", index2_trial, "]")))
#         eval(parse(text = paste0("rm(\"pred_chain_", c, "\")")))
#       }
#
#       pred_trial <- abind::abind(trial_list, along=1)
#       rm(trial_list)
#       eval(parse(text = paste0("pred_trialMCMC_splits[[", s, "]] <- apply(pred_trial, c(1,2), mean, na.rm = TRUE)")))
#       rm(pred_trial)
#
#     }
#
#     pred_trialMCMC_splits <- abind::abind(pred_trialMCMC_splits, along=3)
#     pred_trialMCMC[[t]] <- apply(pred_trialMCMC_splits, c(1,2), mean, na.rm = TRUE)
#     rm(pred_trialMCMC_splits)
#
#   }
#
#   ret <- list()
#
#   ret$pred_subjMCMC <- pred_subjMCMC
#   ret$pred_trialMCMC <- pred_trialMCMC
#
#   return(ret)
#
# }
