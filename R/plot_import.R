#' Plot raw experiment data
#'
#' \code{plot_import} plots a single participant's data from [import_single()],
#' or a single participant (if \code{!is.null(id)}), or all participants' data
#' from [import_multiple()].
#'
#' @param parsed_list [import_single()] or [import_multiple()] output.
#' @param import_single Is the output from [import_single()]?
#' @param id subjID to select if only plots for a single participant are
#' desired. Will also accept a single numeric value i, which will select the
#' ith participant in the output.
#' @param types Types of plot to output, choose from any (or all) of
#' \code{train}, \code{test}, and \code{affect}.
#' @param plt.train List of length <= 2, with the first element a single value
#' or numeric vector of the number of trials to lag in the calculation of
#' cumulative probabilities, and the second element a vector of training types
#' to include.
#' @param plt.test List of length <= 2. The first element specifies the types
#' of test pairs to plot; accepted inputs include "chooseA", "avoidB", "novel",
#' "training", "all", and/or character vector(s) of specific pairs. The second
#' input defines how these options are plotted - either "grouped" or
#' "individual"; note "grouped" plots will not work for custom inputs or if
#' "all" is selected as an option.
#' @param plt.affect List of length <= 2 indicating (1) how many trials to lag
#' (only a single value accepted), and (2) the nouns to plot (can be any of
#' "happy", "confident", "engaged", or "fatigue").
#' @param grp_compare Group to compare on which is found from the participant
#' info. Note that if \code{parsed_list} is split into 2 (i.e., distanced and
#' non-distanced), comparisons will be automatically made on this split.
#' @param grp_names Vector of labels for plot keys for the different groups.
#' An attempt will be made to label these automatically; it is recommended to
#' first leave this list empty to make sure the correct labels are applied.
#' @param recode_na Some grouping variables are \code{NA} in the participant
#' information due to them being asked conditionally.  This option enables these
#' to be recoded as appropriate (e.g., to 0 or \code{FALSE}).
#' @param aff_by_reward Enables affect plots to be compared by whether or not
#' the prior stimulus was rewarded or not. This will override grp_compare (but
#' won't affect other types of plots).
#' @param legend_pos Enables the legend positions to be set manually.
#' @param pal Define a custom colour palette for the plots? Otherwise reverts to
#' defaults.
#' @param font Use a custom font for the plots?
#' @param font_size Base plot font size.
#' @param ... Other arguments, used internally by other functions calling this
#' one.
#'
#' @return Either a single or named \code{list} of \code{ggplot} objects
#'
#' @importFrom magrittr %>%
#' @importFrom rlang := !!
#' @export

plot_import <-
  function(parsed_list,
           import_single = FALSE,
           id = NULL,
           types = c("train", "test", "affect"),
           plt.train = list(),
           plt.test = list(),
           plt.affect = list(),
           grp_compare = NULL,
           grp_names = c(),
           recode_na = NULL,
           aff_by_reward = FALSE,
           legend_pos = "right",
           pal = NULL,
           font = "",
           font_size = 14,
           ...) {

    if (length(grp_compare) > 1) stop("Can only compare on one feature.")
    if (length(grp_names) > 2) warning(
      "Comparisons on features that are not binary are not recommended."
      )

    if (is.null(pal)) {
      pal <- c("#ffc9b5", "#648767", "#b1ddf1", "#95a7ce", "#987284", "#3d5a80")
    }
    else if (!is.null(pal) & length(pal) < 6) {
      message("Need at least 6 colours, reverting to defaults.")
      pal <- c("#ffc9b5", "#648767", "#b1ddf1", "#95a7ce", "#987284", "#3d5a80")
    }

    ## useless assignments to appease R CMD check
    subjID <- exclusion <- group <- type <- value <- choice <- trial_no <-
      trial_no_group <- cuml_acc_mean <- cuml_acc_mean_sub_se <-
      cuml_acc_mean_pl_se <- fatigue_response <- fatigue_slider_start <-
      reward <- question_type <- question_response <- question_slider_start <-
      question_no_group <- cum_resp_l <- cum_resp_l_mean <-
      cum_resp_l_mean_sub_se <- cum_resp_l_mean_pl_se <- test_type <-
      test_type_perf <- count <- mean_prop_correct <- colour_stim <-
      mean_prop_correct_sub_se <- mean_prop_correct_pl_se <- NULL

    l <- list(...)

    if (import_single & is.null(l$test_df)) {
      training <- parsed_list$training
      test <- parsed_list$test
    }
    else if (!import_single & !is.null(id)) {
      if (!is.null(parsed_list$individual_results)) {
        if (is.numeric(id)) id <- names(parsed_list$individual_results)[[id]]
        else id <- paste0("ID", as.character(id))
        training <- parsed_list$individual_results[[id]]$training
        test <- parsed_list$individual_results[[id]]$test
        import_single <- TRUE
      } else {
        if (is.numeric(id)) id <- unique(parsed_list$ppt_info$subjID)[[id]]
        else id <- gsub("#", "", id)
        training <- parsed_list$training %>%
          dplyr::right_join(tibble::as_tibble(id), by = c("subjID" = "value"))
        test <- parsed_list$test %>%
          dplyr::right_join(tibble::as_tibble(id), by = c("subjID" = "value"))
        import_single <- TRUE
      }
    }
    else if (!is.null(parsed_list)) {
      if (length(parsed_list) == 2 & is.null(grp_compare)) {
        ids_vec <- parsed_list[[1]]$ppt_info %>%
          dplyr::bind_rows(parsed_list[[2]]$ppt_info) %>%
          dplyr::select(subjID, exclusion) %>%
          dplyr::filter(exclusion == 0) %>%
          dplyr::select(-exclusion)
        training <- parsed_list[[1]]$training %>%
            dplyr::mutate(group = names(parsed_list)[1]) %>%
            dplyr::bind_rows(
              dplyr::mutate(
                parsed_list[[2]]$training, group = names(parsed_list)[2]
                )
              ) %>%
            dplyr::right_join(ids_vec, by = "subjID") %>%
            dplyr::mutate(
              group = ifelse(
                !is.null(recode_na) & is.na(group), recode_na, group
                )
            )
        test <- parsed_list[[1]]$test %>%
          dplyr::mutate(group = names(parsed_list)[1]) %>%
          dplyr::bind_rows(
            dplyr::mutate(parsed_list[[2]]$test, group = names(parsed_list)[2])
            ) %>%
          dplyr::right_join(ids_vec, by = "subjID") %>%
          dplyr::mutate(
            group = ifelse(!is.null(recode_na) & is.na(group), recode_na, group)
          )
      }
      else if (is.null(grp_compare)) {
        ids_vec <- parsed_list$ppt_info %>%
          dplyr::select(subjID, exclusion) %>%
          dplyr::filter(exclusion == 0) %>%
          dplyr::select(-exclusion)
        training <- parsed_list$training %>%
          dplyr::right_join(ids_vec, by = "subjID")
        test <- parsed_list$test %>%
          dplyr::right_join(ids_vec, by = "subjID")
      }
      else {
        if (length(parsed_list) == 2) {
          ids_vec <- parsed_list[[1]]$ppt_info %>%
            dplyr::bind_rows(parsed_list[[2]]$ppt_info)
        }
        else {
          ids_vec <- parsed_list$ppt_info
        }
        ids_vec <- ids_vec %>%
          dplyr::select(subjID, exclusion, dplyr::any_of(grp_compare)) %>%
          dplyr::filter(exclusion == 0) %>%
          dplyr::select(-exclusion) %>%
          dplyr::rename(group = 2) %>%
          dplyr::mutate(
            group = ifelse(!is.null(recode_na) & is.na(group), recode_na, group)
          )

        if (length(parsed_list) == 2) {
          training <- parsed_list[[1]]$training %>%
            dplyr::bind_rows(parsed_list[[2]]$training) %>%
            dplyr::right_join(ids_vec, by = "subjID")
          test <- parsed_list[[1]]$test %>%
            dplyr::bind_rows(parsed_list[[2]]$test) %>%
            dplyr::right_join(ids_vec, by = "subjID")
        } else {
          training <- parsed_list$training %>%
            dplyr::right_join(ids_vec, by = "subjID")
          test <- parsed_list$test %>%
            dplyr::right_join(ids_vec, by = "subjID")
        }
      }
    }

    ret <- list()

    if (any(types == "train")) {
      pairs <- list("AB", "CD", "EF")
      names(pairs) <- c("12", "34", "56")
      if (!is.null(grp_compare) | length(parsed_list) == 2) {
        if (length(grp_names) == 0) {
          training <- training %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
              type = paste0(pairs[[as.character(type)]], " (", group, ")")
              ) %>%
            dplyr::ungroup()
        }
        else {
          names(grp_names) <- as.character(unique(training$group))
          training <- training %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
              type = paste0(pairs[[as.character(type)]], " (",
                            grp_names[[as.character(group)]], ")")) %>%
            dplyr::mutate(group = grp_names[[as.character(group)]]) %>%
            dplyr::ungroup()
        }
      }
      else {
        training <- training %>%
          dplyr::rowwise() %>%
          dplyr::mutate(type = pairs[[as.character(type)]]) %>%
          dplyr::ungroup()
      }

      if (tryCatch(length(plt.train[[2]]), error = function(e) FALSE)) {
          train_types <- tibble::as_tibble(plt.train[[2]]) %>%
            dplyr::mutate(value = as.character(value))
          train_df <- training %>%
            dplyr::inner_join(train_types, by = c("type" = "value"))
      }
      else if (tryCatch(length(plt.train[[1]]), error = function(e) FALSE)) {
          trial_lags <- plt.train[[1]][is.numeric(plt.train[[1]])]
          train_df <- training %>%
            dplyr::select(-tidyselect::contains("cum_prob")) %>%
            tidyr::drop_na(choice) %>%
            dplyr::arrange(trial_no) %>%
            dplyr::group_by(subjID, type)
          if (!is.null(grp_compare) | length(parsed_list) == 2) {
            train_df <- train_df %>%
              dplyr::group_by(group, .add = TRUE)
          }
          train_df <- train_df %>%
            dplyr::mutate(trial_no_group = dplyr::row_number())
          for (lag in trial_lags) {
            col_name <- rlang::sym(paste0("cuml_accuracy_l", lag))
            train_df <- train_df %>%
              dplyr::mutate(
                !!col_name := runner::runner(
                  x = choice,
                  f = function(x) {sum(x, na.rm = T)/sum(!is.na(x))},
                  k = lag
                  )
                )
          }
          train_df <- train_df %>%
            dplyr::ungroup()
      }
      else {
        train_df <- training
      }
      if (!exists("trial_lags")) {
        trial_lags <- 20
      }

      train_df <- train_df %>%
        dplyr::select(subjID, trial_no_group, type,
                      tidyselect::contains("cuml_accuracy"),
                      dplyr::any_of("group")) %>%
        dplyr::group_by(trial_no_group, type)

      cols <- names(train_df)[startsWith(names(train_df), "cuml_accuracy")]
      tr_plts <- list()
      for (l in seq_along(trial_lags)) {
        n_lag <- trial_lags[l]
        col <- rlang::sym(cols[l])
        plt_name <- paste0("training_lag", n_lag)

        tr_plot_df <- train_df %>%
          dplyr::mutate(cuml_acc_mean = mean(!!col, na.rm = TRUE)) %>%
          dplyr::mutate(cuml_acc_mean_sub_se = cuml_acc_mean - std(!!col)) %>%
          dplyr::mutate(cuml_acc_mean_pl_se = cuml_acc_mean + std(!!col)) %>%
          dplyr::ungroup() %>%
          dplyr::distinct(
            trial_no_group, type, cuml_acc_mean, cuml_acc_mean_sub_se,
            cuml_acc_mean_pl_se
          )

        plt_tr <- tr_plot_df %>%
          ggplot2::ggplot(
            ggplot2::aes(x = trial_no_group, y = cuml_acc_mean,
                         colour = factor(type), fill = factor(type))
            ) +
          ggplot2::geom_point(alpha=0.65) +
          ggplot2::geom_line() +
          ggplot2::scale_x_continuous(breaks = seq(0, 120, 20)) +
          ggplot2::geom_vline(
            xintercept = tryCatch(c(seq(n_lag, 120 - n_lag, n_lag)),
                                  error = function(e) NULL),
            linetype = "dashed", alpha = 0.5
          ) +
          ggplot2::xlab("Trial number") +
          ggplot2::ylab("Cumulative A/C/E choice probability (\u00B1 SE)") +
          ggplot2::scale_color_manual(name = "Trial Type", values = pal) +
          ggplot2::scale_fill_manual(name = "Trial Type",
                                     values = unlist(pal)) +
          cowplot::theme_half_open(
            font_size = font_size,
            font_family = font
          ) +
          ggplot2::theme(legend.position = legend_pos) +
          ggplot2::ggtitle(paste0(n_lag, "-trial lagged"))

        if (!import_single) {
          plt_tr <- plt_tr +
            ggplot2::geom_ribbon(ggplot2::aes(
              ymin = cuml_acc_mean_sub_se, ymax = cuml_acc_mean_pl_se),
              alpha = 0.2
            )
        }

        tr_plts[[plt_name]] <- plt_tr

      }
      ret$training <- tr_plts
    }

    if (any(types == "affect")) {
      if (tryCatch(length(plt.affect[[1]]), error = function(e) FALSE)) {
        a_lag <- plt.affect[[1]]
      }
      else a_lag <- 1
      if (tryCatch(length(plt.affect[[2]]), error = function(e) FALSE)) {
        af_types <- plt.affect[[2]]
      }
      else af_types <- c("happy", "confident", "engaged", "fatigue")

      if (any(af_types == "fatigue")) {
        fatigue_qs <- training %>%
          dplyr::filter(!is.na(fatigue_response)) %>%
          dplyr::mutate(
            dplyr::across(.cols = c(2:3, 8:17), ~ replace(.x, values = NA))
            ) %>%
          dplyr::mutate(question_type = "fatigue") %>%
          dplyr::mutate(question_slider_start = fatigue_slider_start) %>%
          dplyr::mutate(question_response = fatigue_response)
        training <- training %>%
          dplyr::bind_rows(fatigue_qs)
      }

      affect_df <- training %>%
        dplyr::select(subjID, type, choice, reward, trial_no, question_type,
                      question_response, question_slider_start,
                      dplyr::any_of("group")) %>%
        dplyr::arrange(trial_no) %>%
        dplyr::group_by(subjID, question_type) %>%
        dplyr::mutate(
          cum_resp_l =
            runner::runner(
              x = question_response,
              f = function(x) {sum(x, na.rm=T)/sum(!is.na(x))},
              k = a_lag
              )
        ) %>%
        dplyr::mutate(question_no_group = dplyr::row_number()) %>%
        dplyr::group_by(question_no_group, question_type)

      if (aff_by_reward) {
        affect_df <- affect_df %>%
          dplyr::mutate(
            group = ifelse(reward==1, "Rewarded", "Not rewarded")
            ) %>%
          dplyr::group_by(group, .add = TRUE)
      } else if (!is.null(grp_compare) | length(parsed_list) == 2) {
        affect_df <- affect_df %>%
          dplyr::group_by(group, .add = TRUE)
      } else {
        affect_df <- affect_df %>%
          dplyr::mutate(group = question_type)
      }

      af_plts <- list()

      for (n in seq_along(af_types)) {
        noun <- af_types[n]
        df_to_plt <- affect_df %>%
          dplyr::filter(question_type == noun) %>%
          dplyr::mutate(question_type = group) %>%
          dplyr::mutate(cum_resp_l_mean = mean(cum_resp_l, na.rm = TRUE)) %>%
          dplyr::mutate(
            cum_resp_l_mean_sub_se = cum_resp_l_mean - std(cum_resp_l)
            ) %>%
          dplyr::mutate(
            cum_resp_l_mean_pl_se = cum_resp_l_mean + std(cum_resp_l)
            ) %>%
          dplyr::distinct(
            question_type, question_no_group, cum_resp_l_mean,
            cum_resp_l_mean_sub_se, cum_resp_l_mean_pl_se
          )

        af_names <- list("Happiness", "Confidence", "Engagement",
                         "Post-block fatigue")
        names(af_names) <- c("happy", "confident", "engaged", "fatigue")

        if (length(unique(df_to_plt$question_type)) == 1) pal_af <- pal[n]
        else pal_af <- pal

        af_plt <- df_to_plt %>%
          ggplot2::ggplot(
            ggplot2::aes(
              x = question_no_group, y = cum_resp_l_mean,
              color = factor(question_type), fill = factor(question_type)
              )
            ) +
          ggplot2::geom_point(alpha=0.65) +
          ggplot2::geom_line()

        if (noun == "fatigue") {
          af_plt <- af_plt +
            ggplot2::scale_x_continuous(breaks=seq(0,6,1)) +
            ggplot2::geom_vline(
              xintercept=c(seq(1, 5, 1)), linetype="dashed", alpha=0.5
              ) +
            ggplot2::xlab("Block number")

        } else {
          af_plt <- af_plt +
            ggplot2::scale_x_continuous(breaks=seq(0,120,20)) +
            ggplot2::geom_vline(
              xintercept=c(seq(20, 100, 20)), linetype="dashed", alpha=0.5
              ) +
            ggplot2::xlab("Trial number")
        }

        af_plt <- af_plt +
          ggplot2::ylab(paste0(af_names[[noun]], " rating /100")) +
          ggplot2::scale_color_manual(name = NULL, values = pal_af) +
          ggplot2::scale_fill_manual(name = NULL, values = unlist(pal_af)) +
          cowplot::theme_half_open(
            font_size = font_size,
            font_family = font
          ) +
          ggplot2::theme(legend.position = legend_pos) +
          ggplot2::ggtitle(af_names[[noun]])

        if (length(unique(df_to_plt$question_type)) == 1) {
          af_plt <- af_plt +
            ggplot2::guides(colour = "none", fill = "none")
        }

        if (!import_single) {
          af_plt <- af_plt +
            ggplot2::geom_ribbon(ggplot2::aes(
              ymin = cum_resp_l_mean_sub_se, ymax = cum_resp_l_mean_pl_se),
              alpha = 0.2
            )
        }

        af_plts[[noun]] <- af_plt
      }

      ret$affect <- af_plts
    }

    if (any(types == "test")) {

      test_grps <- list(c("AB", "CD", "EF"), c("AC", "AD", "AE", "AF"),
                        c("CB", "DB", "EB", "FB"), c("CE", "CF", "ED", "FD"))
      names(test_grps) <- c("training", "chooseA", "avoidB", "novel")

      all_pairs <- list("AB", "CD", "EF", "AC", "AD", "AE", "AF", "CB", "DB",
                        "EB", "FB", "CE", "CF", "ED", "FD")
      names(all_pairs) <- c("12", "34", "56", "13", "14", "15", "16", "32",
                            "42", "52", "62", "35", "36", "54", "64")

      if (tryCatch(length(plt.test[[2]]), error = function(e) FALSE)) {
        tt_grp <- plt.test[[2]]
      }
      else tt_grp <- "grouped"
      if (tryCatch(length(plt.test[[1]]), error = function(e) FALSE)) {
        test_types_plt <- plt.test[[1]]
      }
      else test_types_plt <- c("chooseA", "avoidB", "novel", "training")

      known_grps <- c("chooseA", "avoidB", "novel", "training")

      if (is.null(parsed_list)) {
        test <- l$test_df
      }

      if (tt_grp == "grouped" &
          (length(setdiff(test_types_plt, known_grps)) == 0)) {
        test <- test %>%
          dplyr::filter(test_type %in% test_types_plt) %>%
          dplyr::mutate(
            test_type = factor(test_type, levels = test_types_plt)
            ) %>%
          dplyr::group_by(subjID, test_type)
      }
      else {
        to_keep <- vector(mode = "character")
        if (any(test_types_plt == "all")) {
          to_keep <- unlist(all_pairs)
        }
        else {
          known_types <- intersect(test_types_plt, known_grps)
          for (t in known_types) {
            to_keep <- c(to_keep, test_grps[[t]])
          }
          to_keep <- c(to_keep, setdiff(test_types_plt, known_grps))
        }

      test <- test %>%
        dplyr::rowwise() %>%
        dplyr::mutate(test_type = all_pairs[[as.character(type)]]) %>%
        dplyr::ungroup() %>%
        dplyr::filter(test_type %in% to_keep) %>%
        dplyr::mutate(test_type = factor(test_type)) %>%
        dplyr::group_by(subjID, test_type)
      }

      if (!is.null(grp_compare) | length(parsed_list) == 2) {
        test <- test %>%
          dplyr::group_by(group, .add = TRUE)
      }

      test_plot_df <- test %>%
        dplyr::mutate(
          test_type_perf = sum(choice, na.rm = TRUE)/length(!is.na(choice))
          ) %>%
        dplyr::select(subjID, choice, test_type, test_type_perf,
                      dplyr::any_of("group")) %>%
        dplyr::ungroup(subjID)

      if (import_single) {
        plot_test <- test_plot_df %>%
          dplyr::ungroup() %>%
          dplyr::mutate(choice = ifelse(is.na(choice), 2, choice)) %>%
          ggplot2::ggplot(ggplot2::aes(
            x = test_type, fill = factor(choice, levels=c(1, 0, 2)))
          ) +
          ggplot2::geom_bar(alpha = 0.7) +
          ggplot2::geom_text(
            stat = "count", family = ifelse(!is.null(font), font, ""),
            ggplot2::aes(
              label = ggplot2::after_stat(count),
              colour = factor(choice, levels = c(1, 0, 2))
            ),
            position = ggplot2::position_stack(vjust=0.5)
          ) +
          ggplot2::xlab("Test type") +
          ggplot2::ylab("Count") +
          ggplot2::scale_fill_manual(
            values = pal, name = NULL,
            labels = c("Correct", "Incorrect", "Timed out")
            ) +
          ggplot2::scale_colour_manual(
            values = c("#000000", "#FFFFFF", "#FF0000"), guide = "none"
            ) +
          cowplot::theme_half_open(
            font_size = font_size,
            font_family = font
          ) +
          ggplot2::ggtitle("Test performance")
      }
      else {
        if (!is.null(grp_compare) | length(parsed_list) == 2) {
          test_plot_df <- test_plot_df %>%
            dplyr::distinct(subjID, test_type, test_type_perf, group)
        }
        else {
          test_plot_df <- test_plot_df %>%
            dplyr::distinct(subjID, test_type, test_type_perf)
        }
        test_plot_df <- test_plot_df %>%
          dplyr::mutate(
            mean_prop_correct = mean(test_type_perf, na.rm = TRUE)
            ) %>%
          dplyr::mutate(
            mean_prop_correct_sub_se = mean_prop_correct - std(test_type_perf)
            ) %>%
          dplyr::mutate(
            mean_prop_correct_pl_se = mean_prop_correct + std(test_type_perf)
            ) %>%
          dplyr::ungroup() %>%
          dplyr::select(-subjID) %>%
          dplyr::distinct()

        if (!is.null(grp_compare) | length(parsed_list) == 2) {
          test_plot_df <- test_plot_df %>%
            dplyr::mutate(colour_stim = factor(group))
        }
        else {
          n_types <- length(unique(test_plot_df$test_type))
          if (length(pal) < n_types) {
            test_plot_df <- test_plot_df %>%
              dplyr::rowwise() %>%
              dplyr::mutate(colour_stim = substr(test_type, 1, 1)) %>%
              dplyr::mutate(
                colour_stim = factor(colour_stim,
                                     levels = c("A", "C", "D", "E", "F"))
                ) %>%
              dplyr::ungroup()
          } else {
              test_plot_df <- test_plot_df %>%
                dplyr::mutate(colour_stim = factor(test_type))
          }
        }

        plot_test <- test_plot_df %>%
          dplyr::distinct(test_type, colour_stim, .keep_all = T) %>%
          ggplot2::ggplot(
            ggplot2::aes(x = test_type, y = mean_prop_correct*100,
                         fill = colour_stim)
            ) +
          ggplot2::geom_bar(stat="identity", alpha = 0.7, color = "black",
                            position = ggplot2::position_dodge()) +
          ggplot2::geom_errorbar(
            ggplot2::aes(
              ymin = mean_prop_correct_sub_se*100,
              ymax = mean_prop_correct_pl_se*100
              ),
            width = .2, position = ggplot2::position_dodge(.9)
          ) +
          ggplot2::xlab("Test type") +
          ggplot2::ylab("% correct (\u00B1 SE)") +
          cowplot::theme_half_open(
            font_size = font_size,
            font_family = font
          ) +
          ggplot2::theme(legend.position = legend_pos) +
          ggplot2::ggtitle("Test performance")

        if (is.null(grp_compare) & length(parsed_list) != 2) {
            plot_test <- plot_test +
              ggplot2::scale_fill_manual(name = NULL, values = pal) +
              ggplot2::guides(fill = "none")
        } else if (length(grp_names) > 0) {
          plot_test <- plot_test +
            ggplot2::scale_fill_manual(name = NULL, values = pal,
                                       labels = grp_names)
        } else {
          plot_test <- plot_test +
              ggplot2::scale_fill_manual(name = NULL, values = pal)
        }

      }

      ret$test$test_perf <- plot_test
    }

  if (length(ret) == 1 & length(ret[[1]]) == 1) ret <- ret[[1]][[1]]
  else if (length(ret) == 1) ret <- ret[[1]]
  return(ret)
}
