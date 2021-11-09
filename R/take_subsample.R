#' Helper function to take a subsample of our parsed list (for demonstration purposes)
#'
#' \code{take_subsample} takes a sample of size n to enable model fitting on a smaller sample.
#'
#' @param parsed_list A list outputted from [import_multiple()].
#' @param n_ppts Sample size to take.
#'
#' @return A named \code{list}.
#'
#' @importFrom magrittr %>%
#' @export

take_subsample <- function(parsed_list, n_ppts) {

  if (is.null(parsed_list[["ppt_info"]])) {
    stop("Could not find a list of participant info to take a sample of IDs. Perhaps the list is split?")
  }

  ids <- sample(unique(parsed_list[["ppt_info"]][["subjID"]]), size = n_ppts)
  subsample <- list()
  elements <- names(parsed_list)

  for (el in elements) {
    subsample[[el]] <- parsed_list[[el]] %>%
      dplyr::filter(subjID %in% ids)
  }

  return(subsample)
}
