#' JATOS text file parsing function (multiple results)
#'
#' \code{import_multiple()} parses a JATOS .txt file containing multiple
#' participant's results to a (list of) R object(s).
#'
#' @param jatos_txt_file JATOS .txt file path.
#' @param separate Returns separate lists for non-distanced and distanced
#' groups?
#' @param exclusion Exclude those meeting exclusion criteria from questionnaire
#' / hBayesDM outputs?
#' @param indiv Return individual-level data? This can then be passed to
#' [plot_import()].
#' and plotted for one individual at a time.
#' @param ret_incomplete Return incomplete datasets? They may require manual parsing.
#' @param ... Other arguments passed to [import_single()].
#'
#' @return \code{list} of \code{tibbles} or a single \code{tibble} if combine = TRUE
#'
#' @importFrom magrittr %>%
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export

import_multiple <- function(jatos_txt_file,
                            separate = TRUE,
                            exclusion = TRUE,
                            indiv = FALSE,
                            ret_incomplete = FALSE,
                            ...) {

  l <- list(...)
  if (is.null(l$hbayesDM)) l$hbayesDM <- FALSE
  if (is.null(l$qstns_gillan)) l$qstns_gillan <- TRUE
  if (is.null(l$prolific)) l$prolific <- TRUE
  if (is.null(l$accuracy)) l$accuracy <- FALSE
  if (is.null(l$task_excl)) l$task_excl <- TRUE
  if (is.null(l$combine)) l$combine <- FALSE
  if (is.null(l$issues)) l$issues <- FALSE
  if (is.null(l$add_sex)) l$add_sex <- FALSE

  message("Reading text file...")

  ## to appease R CMD check
  subjID <- NULL

  all_comp <-
    jsonlite::fromJSON(
      sprintf('[%s]', paste(
        readLines(jatos_txt_file, encoding = "UTF-8", warn = F), collapse = ',')
        )
      )
  res_splits <-
    which(
      lapply(1:length(all_comp),
             function (i) any(
               grepl("Welcome to the experiment.", all_comp[[i]], fixed = T)
               )
             )
      %in% TRUE
    )
    # look along the list of individual component results, note down the list numbers where
    # "Welcome to the experiment." occurs. If all results complete, should be every 5
    # (i.e. 1, 6, 11 etc)

  indiv_res_complete <- list()
  indiv_res_incomplete <- list()

  message("Identifying individual results...")

  for (i in seq_along(res_splits)) {
    end = ifelse(i==length(res_splits),length(all_comp),res_splits[i+1]-1)
    result = all_comp[res_splits[i]:end]
    if (length(result)==5) {
      curr_length = length(indiv_res_complete)
      indiv_res_complete[[curr_length+1]] <- result
    } else {
      curr_length = length(indiv_res_incomplete)
      indiv_res_incomplete[[curr_length+1]] <- result
    }
  }

  if (length(indiv_res_incomplete)>0 & !ret_incomplete) {
    message(
      strwrap(
        "One or more incomplete datasets have been omitted. Set
        incomplete=T to retrieve these data (will require manual parsing).",
        prefix = " ", initial = ""
      )
    )
  }

  if (separate) {
    parsed_results_nd <- list()
    parsed_results_d <- list()

    if (exclusion) {
      excluded_nd <- list()
      included_nd <- list()
      excluded_d <- list()
      included_d <- list()
    }
  } else {
    parsed_results <- list()

    if (exclusion) {
      excluded <- list()
      included <- list()
    }
  }



  message("Parsing results...")

  pb = txtProgressBar(min = 0, max = length(indiv_res_complete),
                      initial = 0, style = 3)

  for (i in seq_along(indiv_res_complete)) {
    setTxtProgressBar(pb, i)

    if (separate) {
      res <- import_single(
        indiv_res_complete[[i]],
        multiple = TRUE,
        accuracy = l$accuracy,
        task_excl = l$task_excl,
        hbayesDM = l$hbayesDM,
        prolific = l$prolific,
        qstns_gillan = l$qstns_gillan,
        combine = l$combine,
        issues = l$issues,
        add_sex = l$add_sex
      )

      if (res$ppt_info$distanced) {
        parsed_results_d[[i]] <- res
        names(parsed_results_d)[i] <- paste0(
          "ID",parsed_results_d[[i]]$ppt_info$subjID)
        if (exclusion) {
          if (parsed_results_d[[i]]$ppt_info$exclusion) {
            excluded_d[[i]] <- parsed_results_d[[i]]
            names(excluded_d)[i] <- paste0(
              "ID", parsed_results_d[[i]]$ppt_info$subjID)
          } else {
            included_d[[i]] <- parsed_results_d[[i]]
            names(included_d)[i] <- paste0(
              "ID",parsed_results_d[[i]]$ppt_info$subjID)
          }
        }
      } else {
        parsed_results_nd[[i]] <- res
        names(parsed_results_nd)[i] <- paste0(
          "ID",parsed_results_nd[[i]]$ppt_info$subjID)
        if (exclusion) {
          if (parsed_results_nd[[i]]$ppt_info$exclusion) {
            excluded_nd[[i]] <- parsed_results_nd[[i]]
            names(excluded_nd)[i] <- paste0(
              "ID", parsed_results_nd[[i]]$ppt_info$subjID)
          } else {
            included_nd[[i]] <- parsed_results_nd[[i]]
            names(included_nd)[i] <- paste0(
              "ID",parsed_results_nd[[i]]$ppt_info$subjID)
          }
        }
      }
    } else {
      parsed_results[[i]] <- import_single(
        indiv_res_complete[[i]],
        multiple = TRUE,
        accuracy = l$accuracy,
        task_excl = l$task_excl,
        hbayesDM = l$hbayesDM,
        prolific = l$prolific,
        qstns_gillan = l$qstns_gillan,
        combine = l$combine,
        issues = l$issues,
        add_sex = l$add_sex
      )
      names(parsed_results)[i] <- paste0(
        "ID",parsed_results[[i]]$ppt_info$subjID)

      if (exclusion) {
        if (parsed_results[[i]]$ppt_info$exclusion) {
          excluded[[i]] <- parsed_results[[i]]
          names(excluded)[i] <- paste0(
            "ID", parsed_results[[i]]$ppt_info$subjID)
        } else {
          included[[i]] <- parsed_results[[i]]
          names(included)[i] <- paste0(
            "ID",parsed_results[[i]]$ppt_info$subjID)
        }
      }
    }
  }

  if (separate) {

    demographics_nd <- list()
    questionnaire_totals_nd <- list()
    training_nd <- list()
    test_nd <- list()

    demographics_d <- list()
    questionnaire_totals_d <- list()
    training_d <- list()
    test_d <- list()

    ret <- list()

    demographics_nd <- lapply(
      seq_along(parsed_results_nd),
      function(i) demographics_nd[[i]] <- parsed_results_nd[[i]]$ppt_info)
    questionnaires_nd <- lapply(
      seq_along(parsed_results_nd),
      function(i)
        questionnaire_totals_nd[[i]] <- parsed_results_nd[[i]]$questionnaires)

    demographics_d <- lapply(
      seq_along(parsed_results_d),
      function(i) demographics_d[[i]] <- parsed_results_d[[i]]$ppt_info)
    questionnaires_d <- lapply(
      seq_along(parsed_results_d),
      function(i)
        questionnaire_totals_d[[i]] <- parsed_results_d[[i]]$questionnaires)

    ret$non_distanced$ppt_info <- dplyr::left_join(
      tibble::as_tibble(
        dplyr::bind_rows(demographics_nd)
        ),
      dplyr::bind_rows(questionnaires_nd),
      by = c("subjID", "sessionID", "studyID")
    )
    ret$distanced$ppt_info <- dplyr::left_join(
      tibble::as_tibble(dplyr::bind_rows(demographics_d)),
      dplyr::bind_rows(questionnaires_d),
      by = c("subjID", "sessionID", "studyID")
    )

    if (exclusion) {
      exclusion_list_nd <- ret$non_distanced$ppt_info %>%
        dplyr::select(subjID, exclusion)
      exclusion_list_d <- ret$distanced$ppt_info %>%
        dplyr::select(subjID, exclusion)
    }

    training_nd <- lapply(
      seq_along(parsed_results_nd),
      function(i) training_nd[[i]] <- parsed_results_nd[[i]]$training
      )
    training_d <- lapply(
      seq_along(parsed_results_d),
      function(i) training_d[[i]] <- parsed_results_d[[i]]$training
      )

    ret$non_distanced$training <- dplyr::bind_rows(training_nd)
    ret$distanced$training <- dplyr::bind_rows(training_d)

    test_nd <- lapply(seq_along(parsed_results_nd),
                      function(i) test_nd[[i]] <- parsed_results_nd[[i]]$test)
    test_d <- lapply(seq_along(parsed_results_d),
                     function(i) test_d[[i]] <- parsed_results_d[[i]]$test)

    ret$non_distanced$test <- dplyr::bind_rows(test_nd)
    ret$distanced$test <- dplyr::bind_rows(test_d)

    if (l$hbayesDM) {
      qlearning_data_nd <- list()
      qlearning_data_d <- list()

      qlearning_data_nd <-
        lapply(
          seq_along(parsed_results_nd),
          function(i)
            qlearning_data_nd[[i]] <- parsed_results_nd[[i]]$qlearning_data
          )
      qlearning_data_d <-
        lapply(
          seq_along(parsed_results_d),
          function(i)
            qlearning_data_d[[i]] <- parsed_results_d[[i]]$qlearning_data
          )

      if (exclusion) {
        qlearning_data_nd <- dplyr::left_join(
          dplyr::bind_rows(qlearning_data_nd), exclusion_list_nd, by="subjID"
        )
        qlearning_data_d <- dplyr::left_join(
          dplyr::bind_rows(qlearning_data_d), exclusion_list_d, by="subjID"
        )

        ret$non_distanced$qlearning_data <- qlearning_data_nd %>%
          dplyr::filter(exclusion==0) %>%
          dplyr::select(-exclusion)
        ret$distanced$qlearning_data <- qlearning_data_d %>%
          dplyr::filter(exclusion==0) %>%
          dplyr::select(-exclusion)
      } else {
        ret$non_distanced$qlearning_data <- dplyr::bind_rows(qlearning_data_nd)
        ret$distanced$qlearning_data <- dplyr::bind_rows(qlearning_data_d)
      }
    }

    if (l$qstns_gillan) {
      gillan_questions_nd <- list()
      gillan_questions_d <- list()

      gillan_questions_nd <-
        lapply(
          seq_along(parsed_results_nd),
          function(i)
            gillan_questions_nd[[i]] <- parsed_results_nd[[i]]$gillan_questions
          )
      gillan_questions_d <-
        lapply(
          seq_along(parsed_results_d),
          function(i)
            gillan_questions_d[[i]] <- parsed_results_d[[i]]$gillan_questions
          )

      if (exclusion) {
        gillan_questions_nd <-
          dplyr::left_join(
            dplyr::bind_rows(
              gillan_questions_nd), exclusion_list_nd, by="subjID"
            )
        gillan_questions_d <-
          dplyr::left_join(
            dplyr::bind_rows(gillan_questions_d), exclusion_list_d, by="subjID"
            )

        gillan_questions_nd <- gillan_questions_nd %>%
          dplyr::filter(exclusion==0) %>%
          dplyr::select(-exclusion)
        gillan_questions_nd[is.na(gillan_questions_nd)] <- 0
        ret$non_distanced$gillan_questions <- gillan_questions_nd

        gillan_questions_d <- gillan_questions_d %>%
          dplyr::filter(exclusion==0) %>%
          dplyr::select(-exclusion)
        gillan_questions_d[is.na(gillan_questions_d)] <-0
        ret$distanced$gillan_questions <- gillan_questions_d

      } else {
        gillan_questions_nd <- dplyr::bind_rows(gillan_questions_nd)
        gillan_questions_nd[is.na(gillan_questions_nd)] <- 0
        ret$non_distanced$gillan_questions <- gillan_questions_nd

        gillan_questions_d <- dplyr::bind_rows(gillan_questions_d)
        gillan_questions_d[is.na(gillan_questions_d)] <- 0
        ret$distanced$gillan_questions <- gillan_questions_d
      }
    }

    if (l$combine) {
      combined_nd <- list()
      combined_d <- list()

      combined_nd <- lapply(
        seq_along(parsed_results_nd),
        function(i) combined_nd[[i]] <- parsed_results_nd[[i]]$full_data
        )
      combined_d <- lapply(
        seq_along(parsed_results_d),
        function(i) combined_d[[i]] <- parsed_results_d[[i]]$full_data
        )

      ret$non_distanced$full_data <- dplyr::bind_rows(combined_nd)
      ret$distanced$full_data <- dplyr::bind_rows(combined_d)
    }

    if (l$issues) {
      issues_nd <- list()
      issues_d <- list()

      issues_nd <- lapply(
        seq_along(parsed_results_nd),
        function(i) issues_nd[[i]] <- parsed_results_nd[[i]]$issues_comments
        )
      issues_d <- lapply(
        seq_along(parsed_results_d),
        function(i) issues_d[[i]] <- parsed_results_d[[i]]$issues_comments
        )

      ret$non_distanced$issues_comments <- tibble::as_tibble(dplyr::bind_rows(issues_nd))
      ret$distanced$issues_comments <- tibble::as_tibble(dplyr::bind_rows(issues_d))
    }

    if (indiv) {
      ret$non_distanced$individual_results <- parsed_results_nd
      ret$distanced$individual_results <- parsed_results_d

      if (exclusion) {
        ret$non_distanced$individual_results$included <-
          included_nd[-which(sapply(included_nd, is.null))]
        ret$non_distanced$individual_results$excluded <-
          tryCatch(
            excluded_nd[-which(sapply(excluded_nd, is.null))],
            error = function(e) message(
              "\nNo excluded non-distanced participants, set exclusion=F.")
            )
        ret$distanced$individual_results$included <-
          included_d[-which(sapply(included_d, is.null))]
        ret$distanced$individual_results$excluded <-
          tryCatch(
            excluded_d[-which(sapply(excluded_d, is.null))],
            error=function(e) message(
              "\nNo excluded distanced participants, set exclusion=F.")
            )
      }
    }

    if (ret_incomplete) {
      ret$incomplete_results <- indiv_res_incomplete
    }

  } else {

    demographics <- list()
    questionnaire_totals <- list()
    training <- list()
    test <- list()

    ret <- list()

    demographics <- lapply(
      seq_along(parsed_results),
      function(i) demographics[[i]] <- parsed_results[[i]]$ppt_info
      )
    questionnaires <- lapply(
      seq_along(parsed_results),
      function(i)
        questionnaire_totals[[i]] <- parsed_results[[i]]$questionnaires
      )
    ret$ppt_info <- dplyr::left_join(
      tibble::as_tibble(dplyr::bind_rows(demographics)),
      dplyr::bind_rows(questionnaires),
      by = c("subjID", "sessionID", "studyID")
    )

    if (exclusion) {
      exclusion_list <- ret$ppt_info %>% dplyr::select(subjID, exclusion)
    }

    training <- lapply(
      seq_along(parsed_results),
      function(i) training[[i]] <- parsed_results[[i]]$training
      )
    ret$training <- dplyr::bind_rows(training)

    test <- lapply(seq_along(parsed_results),
                   function(i) test[[i]] <- parsed_results[[i]]$test)
    ret$test <- dplyr::bind_rows(test)

    if (l$hbayesDM) {
      qlearning_data <- list()
      qlearning_data <- lapply(
        seq_along(parsed_results),
        function(i) qlearning_data[[i]] <- parsed_results[[i]]$qlearning_data
        )
      if (exclusion) {
        qlearning_data <- dplyr::left_join(
          dplyr::bind_rows(qlearning_data), exclusion_list, by="subjID")
        ret$qlearning_data <- qlearning_data %>%
          dplyr::filter(exclusion==0) %>%
          dplyr::select(-exclusion)
      } else {
        ret$qlearning_data <- dplyr::bind_rows(qlearning_data)
      }
    }

    if (l$qstns_gillan) {
      gillan_questions <- list()
      gillan_questions <- lapply(
        seq_along(parsed_results),
        function(i)
          gillan_questions[[i]] <- parsed_results[[i]]$gillan_questions
        )
      if (exclusion) {
        gillan_questions <- dplyr::left_join(
          dplyr::bind_rows(gillan_questions), exclusion_list, by="subjID"
        )
        gillan_questions <- gillan_questions %>%
          dplyr::filter(exclusion==0) %>%
          dplyr::select(-exclusion)

        gillan_questions[is.na(gillan_questions)] <- 0
        ret$gillan_questions <- gillan_questions

      } else {
        gillan_questions <- dplyr::bind_rows(gillan_questions)
        gillan_questions[is.na(gillan_questions)] <- 0
        ret$gillan_questions <- gillan_questions
      }
    }

    if (l$combine) {
      combined <- list()
      combined <- lapply(
        seq_along(parsed_results),
        function(i) combined[[i]] <- parsed_results[[i]]$full_data
        )
      ret$full_data <- dplyr::bind_rows(combined)
    }

    if (l$issues) {
      issues <- list()
      issues <- lapply(
        seq_along(parsed_results),
        function(i) issues[[i]] <- parsed_results[[i]]$issues_comments
        )
      ret$issues_comments <- tibble::as_tibble(dplyr::bind_rows(issues))
    }

    if (indiv) {
      ret$individual_results <- parsed_results
      if (exclusion) {
        ret$individual_results$included <-
          included[-which(sapply(included, is.null))]
        ret$individual_results$excluded <-
          tryCatch(
            excluded[-which(sapply(excluded, is.null))],
            error = function(e) message(
              "\nNo excluded participants, set exclusion=F.")
            )
      }
    }

    if (ret_incomplete) {
      ret$incomplete_results <- indiv_res_incomplete
    }

  }

  return(ret)

}
