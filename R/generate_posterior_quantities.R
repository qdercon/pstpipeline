#' Generate posterior quantities following MCMC
#'
#' \code{generate_posterior_quantities} automates the CmdStanR \code{$generate_quantities} method,
#' which enables the "generated quantities" block of a Stan model to be sampled from after a model is fit.
#' This vastly reduces memory requirements, particularly if per-sample posterior predictions are desired.
#'
#' @param fit_mcmc CmdStanR model environment or a full path to a saved CmdStanR fit (as an .RDS file).
#' @param data_list Raw data that \code{fit} was fit to; either an R object or full path to an .RDS file.
#' @param stan_dir Path to directory where Stan programmes are found.
#' @param out_dir Path to output directory.
#' @param save_model_as Name to give model; defaults to model name from the fit metadata.
#' @param return Return something other than the fit environment? Options are "paths" (for .csv file paths), or
#' "draws_list" which returns only the posterior predictions as a [posterior::draws_list].
#' @param par_chains Maximum number of chains to compute in parallel; defaults to \code{options(mc.cores)} if
#' this has been set, or 4 if not.

generate_posterior_quantities <-
  function(fit_mcmc, data_list, stan_dir = "stan_files/", out_dir = "cmdstan_output/",
           save_model_as = "", return = "", par_chains = getOption("mc.cores", 4)) {

  if (!exists(as.character(fit_mcmc))) fit_mcmc <- readRDS(fit_mcmc)
  if (!exists(as.character(data_list))) data_cmdstan <- readRDS(fit_mcmc)

  if (grepl("_ppc", fit_mcmc$metadata()$model_name)) {
    warning("Model name contains '_ppc' - are you sure predictions don't already exist for this model?")
  }
  stan_file <- paste0(stan_dir, gsub("_model|_ppc", "", fit_mcmc$metadata()$model_name), "_gq", ".stan")
  stan_model <- cmdstan_model(stan_file)

  fit_gq <- stan_prog$generate_quantities(
    fitted_params = fit_mcmc,
    data = data_cmdstan,
    output_dir = out_dir,
    parallel_chains = par_chains
  )

  ## rename csv output files
  outnames <- fit_mcmc$output_files()
  if (save_model_as == "") save_model_as <- fit_mcmc$metadata()[["model_name"]]
  csv_files <- vector(mode = "character", length = length(outnames))

  for (o in seq_along(outnames)) {
    chain_no <- strsplit(outnames[o], "-")[[1]][3]
    csv_files[o] <-
      paste0(out_dir, save_model_as,
             paste0("_gq_", fit_mcmc$metadata()$iter_sampling * fit_mcmc$metadata()$chains,
                    "chain_", chain_no, ".csv"))
    file.rename(from = outnames[o], to = csv_files[o])
  }

  if (return == "paths") return(csv_files)
  else if (return == "draws_list") return(fit_gq$draws(variables = "y_pred", format = "list"))
  else return(fit_gq)
}
