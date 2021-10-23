#' Generate posterior quantities following MCMC
#'
#' \code{generate_posterior_quantities} automates the CmdStanR \code{$generate_quantities} method,
#' which enables the "generated quantities" block of a Stan model to be sampled from after a model is fit.
#' This vastly reduces memory requirements, particularly if per-sample posterior predictions are desired.
#'
#' @param fit_mcmc CmdStanR model environment or a full path to a saved CmdStanR fit (as an .RDS file).
#' @param data_list Raw data that \code{fit} was fit to; either an R object or full path to an .RDS file.
#' @param out_dir Path to output directory.
#' @param save_model_as Name to give model; defaults to model name from the fit metadata.
#' @param return Return something other than the fit environment? Options are "paths" (for .csv file paths), or
#' "draws_list" which returns only the posterior predictions as a [posterior::draws_list].
#' @param par_chains Maximum number of chains to compute in parallel; defaults to \code{options(mc.cores)} if
#' this has been set, or 4 if not.

generate_posterior_quantities <-
  function(fit_mcmc, data_list, out_dir = "outputs/cmdstan/", save_model_as = "",
           return_type = "", par_chains = getOption("mc.cores", 4)) {

  if (grepl("_ppc", fit_mcmc$metadata()$model_name)) {
    warning("Are you sure predictions don't already exist for this model?")
  }

  if (grepl("test", fit_mcmc$metadata()$model_name)) train_test <- "test"
  else train_test <- "training"

  if (grepl("gainloss", fit_mcmc$metadata()$model_name)) alphas <- "2a"
  else alphas <- "1a"

  stan_model <- cmdstanr::cmdstan_model(
    system.file(
      paste(
        "extdata/stan_files/pst", ifelse(alphas == "2a", "gainloss_Q", "Q"), train_test, "gq.stan", sep = "_"
        ),
      package = "pstpipeline"
      ),
    )

  fit_gq <- stan_model$generate_quantities(
    fitted_params = fit_mcmc,
    data = data_list,
    output_dir = out_dir,
    parallel_chains = par_chains
  )

  ## rename csv output files
  outnames <- fit_gq$output_files()
  if (save_model_as == "") {
    save_model_as <- paste(
      fit$metadata()[["model_name"]], model, ifelse(test, "test", "training"), "gq", sep = "_"
    )
  }
  csv_files <- vector(mode = "character", length = length(outnames))

  for (o in seq_along(outnames)) {
    chain_no <- strsplit(outnames[o], "-")[[1]][3]
    csv_files[o] <-
      paste0(out_dir, save_model_as,
             paste0("_gq_", fit_gq$metadata()$iter_sampling * fit_gq$metadata()$chains,
                    "chain_", chain_no, ".csv"))
    file.rename(from = outnames[o], to = csv_files[o])
  }

  if (return_type == "paths") return(csv_files)
  else if (return_type == "draws_list") return(fit_gq$draws(variables = "y_pred", format = "list"))
  else return(fit_gq)
}