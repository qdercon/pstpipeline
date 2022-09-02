#' JATOS text file parsing function (single result)
#'
#' \code{import_single()} parses a single JATOS .txt file to a (list of) R
#' object(s).
#'
#' @param jatos_txt_file JATOS text file path
#' @param accuracy Include final block AB accuracy < 0.6 as part of exclusion
#' criteria?
#' @param task_excl Use digit span of 0 or > 95% preference for one button as
#' exclusion criteria?
#' @param hbayesDM Return additional df with training data formatted for
#' hBayesDM?
#' @param qstns_gillan Return formatted subset of questionnaire questions?
#' @param prolific Use prolific IDs as subjID?
#' @param combine Combine everything into a single data.frame?
#' @param issues Return feedback from end of experiment?
#' @param incomplete Is the text file of interest incomplete?
#' @param add_sex Add sex to participant info? Imputes (assumed) birth-assigned
#' sex for non-binary individuals from Prolific export - workaround to enable it
#' to be more reasonably controlled for in models (given very low numbers of
#' non-binary individuals).
#' @param ... Other arguments (used by [import_multiple()]).
#'
#' @returns \code{list} of \code{tbl_df}, or a single \code{tbl_df} if
#' \code{combine == TRUE}.
#'
#' @export

import_single <- function(jatos_txt_file,
                          accuracy = FALSE,
                          task_excl = TRUE,
                          hbayesDM = FALSE,
                          qstns_gillan = TRUE,
                          prolific = TRUE,
                          combine = FALSE,
                          issues = FALSE,
                          incomplete = FALSE,
                          add_sex = FALSE,
                          ...) {

  l <- list(...)
  if (is.null(l$multiple)) l$multiple <- FALSE

  if (l$multiple | incomplete) {
    all_comp = jatos_txt_file
  } else {
    all_comp <- jsonlite::fromJSON(
      sprintf('[%s]', paste(
        readLines(jatos_txt_file, encoding="UTF-8", warn=F), collapse = ',')
        )
      )
    if (length(all_comp)!=5) {
      stop(
        "Input is a non-standard length - check for missing tasks in text file."
        )
    }
  }

  if (prolific) {
    pro_id <- list(tryCatch(all_comp[[3]]$prolific_id[1],
                            error = function(e) return(NA)),
                   tryCatch(all_comp[[4]]$prolific_id[1],
                            error = function(e) return(NA)),
                   tryCatch(all_comp[[5]]$prolific_id[1],
                            error = function(e) return(NA)))
    pro_id <- pro_id[!is.na(pro_id)]

    stu_id <- list(tryCatch(all_comp[[3]]$study_id[1],
                            error = function(e) return(NA)),
                   tryCatch(all_comp[[4]]$study_id[1],
                            error = function(e) return(NA)),
                   tryCatch(all_comp[[5]]$study_id[1],
                            error = function(e) return(NA)))
    stu_id <- stu_id[!is.na(stu_id)]

    ses_id <- list(tryCatch(all_comp[[3]]$session_id[1],
                            error = function(e) return(NA)),
                   tryCatch(all_comp[[4]]$session_id[1],
                            error = function(e) return(NA)),
                   tryCatch(all_comp[[5]]$session_id[1],
                            error = function(e) return(NA)))
    ses_id <- ses_id[!is.na(ses_id)]

    if (!all(sapply(pro_id, FUN = identical, all_comp[[1]]$prolific_id[1]))) {
      warning(paste0("Prolific IDs inconsistent, recording as NA: ",
                     paste(pro_id, sep=", ")))
      subjID=NA
    } else {
      subjID=all_comp[[1]]$prolific_id[1]
    }
    if (!all(sapply(stu_id, FUN = identical, all_comp[[1]]$study_id[1]))) {
      warning(paste0("Study IDs inconsistent, recording as NA: ",
                     paste(stu_id, sep=", ")))
      studyID=NA
    } else {
      studyID=all_comp[[1]]$study_id[1]
    }
    if (!all(sapply(ses_id, FUN = identical, all_comp[[1]]$session_id[1]))) {
      warning(paste0("Session IDs inconsistent, recording as NA: ",
                     paste(ses_id, sep=", ")))
      sessionID=NA
    } else {
      sessionID=all_comp[[1]]$session_id[1]
    }
  }
  else {
    subjID=sample(1:100000, 1)
    studyID=NA
    sessionID=NA
  }

  if (!l$multiple) message(paste0("Importing data for subjID #", subjID))

  distanced = data.frame(all_comp[[2]]=="B")
  names(distanced) <- "distanced"

  ret <- list()

  ## useless assignments to appease R CMD check
  age <- question <- questionnaire <- question_no <- catch_question_pass <-
    ZDS_total <- OCIR_total <- LSAS_fear_total <- LSAS_avoidance_total <-
    EAT_total <- BIS_total <- STAI_total <- AES_total <- DARS_hobbies <-
    DARS_food_drink <- DARS_social <- DARS_sensory <- DARS_total <-
    MFIS_physical <- MFIS_cognitive <- MFIS_psychosocial <- MFIS_total <-
    BPQ_total <- SPQ_cognitive_perceptual <- SPQ_interpersonal <-
    SPQ_disorganised <- SPQ_total <- gender <- test_part <- trial_index <-
    trial_type <- correct <- timeout <- stimulus <- glyph_seq <- type <-
    trial_no <- choice <- reward <- points <- trial_block <- trial_no_group <-
    key_press <- correct_response <- cuml_accuracy_l20 <- slider_start <-
    response <- question_type <- question_slider_start <- question_rt <-
    question_response <- trial_time <- time_elapsed <- fatigue_rt <-
    fatigue_slider_start <- fatigue_response <- participant_id <- Sex <-
    sex_prolific <- test_type <- test_trial_no_group <- test_trial_no <-
    cuml_accuracy_test <- gillan_name <- score <- enjoyment <- concentration <-
    NULL


  ## demographics

  demographics <- all_comp[[1]] |>
    dplyr::select(age:length(all_comp[[1]]), -question,
                  -tidyselect::contains("_id", ignore.case = F)) |>
    dplyr::summarise(dplyr::across(.fns=~max(.x, na.rm=T)))

  if (!is.null(demographics$neurological_disorder)) {
    if (demographics$neurological_disorder=="None") {
      demographics$neurological = 0
      demographics$neurological_disorder = NA
    } else {
      demographics$neurological_disorder <-
        paste(unlist(jsonlite::parse_json(
          demographics$neurological_disorder)
          ), collapse=', ')
    }
  }

  if (!is.null(demographics$other_psych_neurdev_condition)) {
    if (demographics$other_psych_neurdev_condition=="None") {
      demographics$psych_neurdev = 0
      demographics$psych_neurdev_condition = NULL
      demographics$other_psych_neurdev_condition = NULL
    }
  }

  if (!is.null(demographics$psych_neurdev_condition)) {
      demographics$psych_neurdev_condition <-
        paste(unlist(jsonlite::parse_json(
          demographics$psych_neurdev_condition)
          ), collapse=', ')
  }

  if (any(demographics$antidepressant)) {
    demographics$antidepressant_type <-
      paste(unlist(jsonlite::parse_json(
        demographics$antidepressant_type)
        ), collapse=', ')
  }

  if (any(demographics$other_medication)) {
    demographics$other_medication_type <-
      paste(unlist(jsonlite::parse_json(
        demographics$other_medication_type)
        ), collapse=', ')
  }

  if (any(!demographics$long_covid == "No") &
      any(demographics$long_covid_symptoms != "None")) {
    demographics$long_covid_symptoms <-
      paste(unlist(jsonlite::parse_json(
        demographics$long_covid_symptoms)
        ), collapse=', ')
  }

  digit_span = data.frame(all_comp[[4]]$final_digit_span) |> tidyr::drop_na()
  names(digit_span) <- "digit_span"

  questionnaires <- all_comp[[5]]

  catch_questions <- questionnaires |>
    dplyr::filter(questionnaire == "catch_questions") |>
    dplyr::select(question_no, catch_question_pass) |>
    dplyr::arrange(question_no) |>
    tidyr::pivot_wider(names_from = question_no,
                       values_from = catch_question_pass,
                       names_prefix = "catch_question_")

  questionnaires_tr <- questionnaires |>
    dplyr::select(
      ZDS_total, OCIR_total, LSAS_fear_total, LSAS_avoidance_total, EAT_total,
      BIS_total, STAI_total, AES_total, tidyselect::any_of("AUDIT_score"),
      DARS_hobbies, DARS_food_drink, DARS_social, DARS_sensory, DARS_total,
      MFIS_physical, MFIS_cognitive, MFIS_psychosocial, MFIS_total, BPQ_total,
      SPQ_cognitive_perceptual, SPQ_interpersonal, SPQ_disorganised,
      SPQ_total) |>
    dplyr::summarise(dplyr::across(.fns=~max(.x, na.rm=T)))

  ppt_info <- cbind(subjID, sessionID, studyID, distanced, digit_span,
                    catch_questions, demographics)

  if (add_sex & prolific) id_gender <- ppt_info |>
    dplyr::select(subjID, gender)

  questionnaires_tr <- cbind(subjID, sessionID, studyID, questionnaires_tr)


  ## MAIN TASK

  main_task <- cbind(subjID, all_comp[[3]])
  total_points = as.integer(max(main_task$total_points, na.rm=T))

  # useful indexes

  end_block_ind <- main_task |>
    dplyr::filter(test_part == "fatigue_question") |>
    dplyr::select(trial_index)

  final_training_ind <- main_task |>
    dplyr::filter(test_part == "training_complete") |>
    dplyr::select(trial_index)

  # remove uninformative variables

  main_task <- main_task |>
    dplyr::filter(trial_type != "html-button-response",
                  trial_type != "video-button-response",
                  trial_type != "instructions",
                  test_part != "fixation",
                  test_part != "fixation",
                  test_part != "fixation_test",
                  test_part != "feedback",
                  test_part != "congrats",
                  test_part != "training_complete",
                  test_part != "test_quiz",
                  !is.na(test_part)) |>
    dplyr::select(-tidyselect::contains("_id", ignore.case=F))

  training_blocks <- main_task |>
    dplyr::filter(trial_index < final_training_ind[1,1])
  test_block <- main_task |>
    dplyr::filter(trial_index > final_training_ind[1,1])

  ## TRAINING BLOCKS

  # add block numbers

  training_blocks <- training_blocks |>
    dplyr::mutate(trial_block = ifelse(trial_index <= end_block_ind[1,1], 1,
         ifelse(trial_index <= end_block_ind[2,1], 2,
                ifelse(trial_index <= end_block_ind[3,1], 3,
                       ifelse(trial_index <= end_block_ind[4,1], 4,
                              ifelse(trial_index <= end_block_ind[5,1], 5,
                                     ifelse(trial_index <= end_block_ind[6,1],
                                            6, NA)))))))
  # training trials

  training_trials <- training_blocks |>
    dplyr::filter(!grepl("question", test_part)) |>
    dplyr::mutate(type =
             ifelse(grepl("A", test_part), 12,
                    ifelse(grepl("C", test_part), 34,
                           ifelse(grepl("E", test_part), 56, NA)))) |>
    dplyr::mutate(choice =
             ifelse(
               ((test_part=="AB" | test_part=="BA" | test_part=="CD" |
                   test_part=="DC" | test_part=="EF" |
                   test_part=="FE") & correct==T) |
                 # chose choice A/C/E on a standard trial
               ((test_part=="AB_rev" | test_part=="BA_rev" |
                   test_part=="CD_rev" | test_part=="DC_rev" |
                   test_part=="EF_rev" | test_part=="FE_rev") &
                  correct==F & !timeout
                ), 1,
                # chose choice A/C/E on a flipped trial
              ifelse(!timeout, 0, NA))) |>
    dplyr::mutate(reward = ifelse(timeout, NA, as.numeric(correct))) |>
      # correct==F on timeout trials, so need extra check
    dplyr::mutate(stimulus = gsub("<.*?>", "", stimulus)) |>
    dplyr::mutate(glyph_seq = paste(glyph_seq[[1]], collapse='')) |>
    dplyr::group_by(type) |>
    dplyr::arrange(trial_no) |>
    dplyr::mutate(trial_no_group = dplyr::row_number()) |>
    dplyr::mutate(
      cuml_accuracy_l20 =
        runner::runner(x=choice,
                       f=function(x) {sum(x, na.rm=T)/sum(!is.na(x))},
                       k=20)
      ) |>
    dplyr::ungroup() |>
    dplyr::select(subjID, type, choice, reward, trial_block, trial_no,
                  trial_no_group, glyph_seq, stimulus, test_part, key_press, rt,
                  correct_response, correct, timeout, points, cuml_accuracy_l20)


  # training questions

  training_questions <- training_blocks |>
    dplyr::filter(test_part=="question") |>
    dplyr::arrange(trial_index) |>
    dplyr::mutate(trial_no = dplyr::row_number()) |>
    dplyr::rename(question_rt = rt, question_slider_start = slider_start,
                  question_response = response, trial_time = time_elapsed) |>
    dplyr::mutate(stimulus = gsub("<.*?>", "", stimulus)) |>
    dplyr::mutate(trial_time = (trial_time - min(trial_time)) / 60000) |>
    dplyr::mutate(
      question_type = ifelse(grepl("happy", stimulus), "happy",
                             ifelse(grepl("engaged", stimulus), "engaged",
                                    ifelse(grepl("confident", stimulus),
                                           "confident", "fatigue")))) |>
    dplyr::select(subjID, trial_block, trial_no, trial_time, question_type,
                  question_slider_start, question_rt, question_response)

  fatigue_questions <- training_blocks |>
    dplyr::filter(test_part=="fatigue_question") |>
    dplyr::mutate(trial_no = trial_block*60) |>
    dplyr::rename(fatigue_rt = rt, fatigue_slider_start = slider_start,
                  fatigue_response = response) |>
    dplyr::select(subjID, trial_block, trial_no, fatigue_rt,
                  fatigue_slider_start, fatigue_response)

  training <- training_trials |>
    dplyr::left_join(training_questions,
                     by = c("subjID", "trial_block", "trial_no")) |>
    dplyr::left_join(fatigue_questions,
                     by = c("subjID", "trial_block", "trial_no"))

  final_block_AB <- training |>
    dplyr::filter(trial_block==6 & trial_no_group==120 & type==12) |>
    dplyr::select(cuml_accuracy_l20)

  final_block_CD <- training |>
    dplyr::filter(trial_block==6 & trial_no_group==120 & type==34) |>
    dplyr::select(cuml_accuracy_l20)

  final_block_EF <- training |>
    dplyr::filter(trial_block==6 & trial_no_group==120 & type==56) |>
    dplyr::select(cuml_accuracy_l20)

  final_block_AB = final_block_AB[[1]]
  final_block_CD = final_block_CD[[1]]
  final_block_EF = final_block_EF[[1]]

  num_f = sum(training$key_press==70, na.rm=T)
  num_j = sum(training$key_press==74, na.rm=T)
  keypress_percent = max((num_f/(num_f+num_j)*100), (num_j/(num_f+num_j)*100))

  if (!incomplete) {
    if (ppt_info$english=="A1/A2" | ppt_info$english=="B1" |
        ppt_info$neurological |
       !ppt_info$catch_question_1 | !ppt_info$catch_question_2 |
       (!ppt_info$catch_question_3 & !ppt_info$catch_question_4) |
       (final_block_AB < 0.6 & accuracy) |
       (digit_span==0 & task_excl) | (keypress_percent>95 & task_excl)) {
          ppt_info <- cbind(subjID, sessionID, studyID, distanced)
          ppt_info$exclusion=1
    }
    else {
      ppt_info <- cbind(subjID, sessionID, studyID, distanced)
      ppt_info$exclusion=0
    }

    time_taken = sum(c(max(all_comp[[1]]$time_elapsed),
                     max(all_comp[[3]]$time_elapsed),
                     max(all_comp[[4]]$time_elapsed),
                     max(all_comp[[5]]$time_elapsed)))/60000
  }
  else {
    ppt_info <- cbind(subjID, sessionID, studyID, distanced)
    ppt_info$exclusion=NA
  }

  ppt_info$final_block_AB = final_block_AB
  ppt_info$final_block_CD = final_block_CD
  ppt_info$final_block_EF = final_block_EF
  ppt_info$total_points = total_points

  if (!incomplete) {
    ppt_info$total_time_taken = time_taken
  } else {
    ppt_info$total_time_taken = NA
  }

  ppt_info$keypress_percent = keypress_percent
  ppt_info$mean_rt = mean(training$rt, na.rm=T)

  if (add_sex & prolific) {
    id_sex <- exprtd_dmgrphcs |>
      dplyr::select(participant_id, Sex) |>
      dplyr::rename(subjID = participant_id, sex_prolific = Sex) |>
      dplyr::right_join(id_gender, by = "subjID") |>
      dplyr::mutate(sex = ifelse(gender == "Non-binary", sex_prolific, gender))
    sex <- id_sex[["sex"]]
    demographics <- cbind(sex, demographics)
  }


  if (l$multiple) {
    ret$ppt_info <- tibble::as_tibble(
      cbind(ppt_info, digit_span, catch_questions, demographics)
      )
    ret$questionnaires <- tibble::as_tibble(questionnaires_tr)
  } else {
    ppt_info <- cbind(ppt_info, digit_span, catch_questions, demographics)
    ret$ppt_info <- tibble::as_tibble(
      dplyr::left_join(
        ppt_info, questionnaires_tr, by=c("subjID", "sessionID", "studyID")
        )
      )
  }

  ret$training <- training

  ## TEST BLOCK
  chars <-
    unique(
      test_block$test_part
      )[!unique(test_block$test_part) %in% "question_test"]
  make_type_dict <- function(test_types) {
    vals <- list(1, 2, 3, 4, 5, 6)
    names(vals) <- c("A", "C", "E", "F", "D", "B")
    types <- list(1, 3, 5, 6, 4, 2)
    names(types) <- seq(1,6)
    type <- list()

    for (t in seq_along(test_types)) {
      test <- test_types[t]
      type[[t]] <- paste0(types[[as.character(min(vals[[substr(test, 1, 1)]],
                                                vals[[substr(test, 2, 2)]]))]],
                          types[[as.character(max(vals[[substr(test, 1, 1)]],
                                                  vals[[substr(test, 2, 2)]]))]]
                          )
    }
    names(type) <- test_types
    return(type)
  }
  type_key <- make_type_dict(chars)
  get_type_val <- function(test_part, key) {
    return(key[[test_part]])
  }

  test_trials <- test_block |>
    dplyr::filter(!grepl("question", test_part)) |>
    dplyr::mutate(test_type = ifelse(grepl("(AB|BA)|(CD|DC)|(EF|FE)", test_part),
                                     "training",
                                  ifelse(grepl("A", test_part), "chooseA",
                                     ifelse(grepl("B", test_part), "avoidB",
                                            "novel")))) |>
    dplyr::rowwise() |>
    dplyr::mutate(type = as.numeric(get_type_val(test_part, type_key))) |>
    dplyr::mutate(choice = ifelse(!timeout, as.numeric(correct), NA)) |>
    dplyr::group_by(test_type) |>
    dplyr::mutate(test_trial_no_group = dplyr::row_number()) |>
    dplyr::mutate(cuml_accuracy_test = cumsum(correct)/test_trial_no_group) |>
    dplyr::ungroup() |>
    dplyr::mutate(stimulus = gsub("<.*?>", "", stimulus)) |>
    dplyr::mutate(glyph_seq = paste(glyph_seq[[1]], collapse='')) |>
    dplyr::select(subjID, type, choice, test_type, test_trial_no,
                  test_trial_no_group, rt, glyph_seq, stimulus, test_part,
                  key_press, correct_response, correct, timeout,
                  cuml_accuracy_test)

  test_questions <- test_block |>
    dplyr::filter(test_part == "question_test") |>
    dplyr::arrange(trial_index) |>
    dplyr::mutate(test_trial_no = dplyr::row_number()) |>
    dplyr::rename(question_rt = rt, question_slider_start = slider_start,
                  question_response = response, trial_time = time_elapsed) |>
    dplyr::mutate(trial_time = (trial_time - min(trial_time)) / 60000) |>
    dplyr::mutate(stimulus = gsub("<.*?>", "", stimulus)) |>
    dplyr::mutate(question_type =
                    ifelse(grepl("happy", stimulus), "happy",
                           ifelse(grepl("engaged", stimulus), "engaged",
                                  ifelse(grepl("confident", stimulus),
                                         "confident", "fatigue")))) |>
    dplyr::select(subjID, test_trial_no, trial_time, question_type,
                  question_slider_start, question_rt, question_response)

  test <- test_trials |>
    dplyr::left_join(test_questions, by=c("subjID", "test_trial_no")) |>
    dplyr::rename(trial_no = test_trial_no)

  ret$test <- test

  if (hbayesDM) {
    qlearning_data <- training[1:4]
    ret$qlearning_data <- qlearning_data |>
      tidyr::drop_na()
  }

  if (qstns_gillan) {

    template <- as.data.frame(
    matrix(nrow=0, ncol=210,
     dimnames =
       list(NULL,
            c("subjID",
              sapply(1:43, function(i) paste0("SCZ_",i)),
              sapply(1:18, function(i) paste0("OCI_",i)),
              sapply(1:26, function(i) paste0("EAT_",i)),
              sapply(1:18, function(i) paste0("AES_",i)),
              sapply(1:10, function(i) paste0("AUDIT_",i)),
              sapply(1:20, function(i) paste0("SDS_",i)),
              sapply(1:20, function(i) paste0("STAI_",i)),
              sapply(1:30, function(i) paste0("BIS_",i)),
              sapply(1:24, function(i) paste0("LSAS_",i))
              )
            )
     )
    )

    if (prolific) {
      template$subjID <- as.character(template$subjID)
    }

    qstns <- questionnaires |>
      dplyr::filter(!is.na(gillan_name)) |>
      dplyr::select(gillan_name, score) |>
      tidyr::pivot_wider(
        names_from = gillan_name, values_from = score, values_fn = sum
        )

    qstns <- cbind(subjID, qstns)
    ret$gillan_questions <- dplyr::bind_rows(template, qstns)

  }

  if (combine) {
    training <- training |> dplyr::mutate(task_part = "training")
    test <- test |> dplyr::mutate(task_part = "test")
    full_task <- training |> dplyr::bind_rows(test)
    ret$full_data <- ppt_info |> dplyr::left_join(full_task, by="subjID")
  }

  if (issues) {
    issues_comments <- questionnaires |>
      dplyr::select(enjoyment, concentration,
                    tidyselect::contains("technical"),
                    tidyselect::contains("layout"),
                    tidyselect::contains("other")) |>
      tidyr::drop_na()
    ret$issues_comments <- tibble::as_tibble(issues_comments)
  }
  return(ret)
}
