# pstpipeline

[![R-CMD-check](https://github.com/qdercon/pstpipeline/actions/workflows/main.yml/badge.svg)](https://github.com/qdercon/pstpipeline/actions/workflows/main.yml)
[![Docker](https://img.shields.io/badge/dockerhub-image-important.svg?logo=Docker)](https://hub.docker.com/repository/docker/qdercon/pstpipeline/general#)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/qdercon/pstpipeline/blob/main/LICENSE)

**p**robabilistic **s**election **t**ask **pipeline**

## an R package to clean, analyse, and present data from a large online learning study

**Full methods and results from this study can be found in the paper** <small>(† = equal contribution)</small>:

> Dercon†, Q., Mehrhof†, S. Z., Sandhu, T. R., Hitchcock, C., Lawson, R. P., Pizzagalli, D. A., Dalgleish, T., & Nord, C. L. (2023). A core component of psychological therapy causes adaptive changes in computational learning mechanisms. *Psychological Medicine*, 1–11. https://doi.org/10.1017/S0033291723001587

All analyses in the paper can be visually inspected (and, in theory, re-run) via the following Jupyter notebooks, including model fitting and checks:

1. Data cleaning, transdiagnostic psychiatric symptom factor derivation, and plotting of behavioural data: [```data_cleaning_factor_derivation.ipynb```](https://github.com/qdercon/pstpipeline/blob/main/notebooks/data_cleaning_factor_derivation.ipynb).
2. Fitting of all computational models, plus model checks, plots of posterior predictions, and parameter recovery: [```model_fitting_mcmc.ipynb```](https://github.com/qdercon/pstpipeline/blob/main/notebooks/model_fitting_mcmc.ipynb). (An [additional notebook](https://github.com/qdercon/pstpipeline/blob/main/notebooks/model_fitting_vb.ipynb) with models fitted via approximate inference is also provided, which can be far more easily re-run.)
3. Outcome analyses including associations between model parameters and transdiagnostic symptom factors and the distancing intervention: [```main_results.ipynb```](https://github.com/qdercon/pstpipeline/blob/main/notebooks/main_results.ipynb).

### Why an R package?

This package, which accompanies the paper, is not meant to be a brand new toolkit &mdash; better, more flexible packages are available for analysing computational psychiatry experiments. Indeed, some of the ```R``` and ```Stan``` code is inspired by and/or modified from other ```R``` packages, in particular [```hBayesDM```](https://github.com/CCS-Lab/hBayesDM) [[1](#references)] and [```rstanarm```](https://mc-stan.org/rstanarm/).

Instead, its main aims are as follows:

1.  To make it easier for our specific analyses to be replicated by others without lengthy scripts and function definitions &mdash; the package loads all necessary dependencies and custom functions (see below) in the background.
2.  To demonstrate a complete pre- and post-processing pipeline for a common learning task, which (hopefully) shows that such workflows are a) not overwhelmingly difficult to adopt, and b) can elicit valuable mechanistic insights.
3.  To do the above in a high-level manner, while still giving the user control over key aspects - most functionality of the package can be achieved with single-line function calls.

**Update (06/23)**: The package has been extensively updated to include a number of extended $Q$-learning models which include trial-by-trial affect ratings (participants rated either their subjective happiness, confidence, or engagement after each trial). These models are versions of a computational model of subjective happiness derived by Rutledge *et al.* (2014) [[2](#references)], with additional parameters to capture time-dependent affective declines (recently termed "mood drift over time" [[3](#references)]). Most functions described below have been extensively updated to accomodate these models, all of which have examples (see documentation). Full details on the modelling and rationale can be found in the accompanying Jupyter notebook [```affect_model_vb.ipynb```](https://github.com/qdercon/pstpipeline/blob/main/notebooks/affect_model_vb.ipynb). Note that MCMC is not possible for these models due to computational complexity, so approximate inference is used throughout.

## Using the package

### Local R installation

To install the R package and all dependencies directly, run the following:

```R
# install.packages("remotes")
remotes::install_github("qdercon/pstpipeline")
```

The majority of the code in the notebooks is in R, so if you wish to run things this way I would recommend taking a look at the notebooks, and then copy/paste or write your own function calls as appropriate. All user functions (listed below) are fully documented; this documentation can be easily accessed via the ```?``` prefix in R/RStudio. Though written primarily for our specific data/analyses, the functions are written to be relatively flexible, and aspects can be easily modified (e.g., to add new models).

### Docker container

A Docker container containing all package dependencies is also provided - see this [README](https://github.com/qdercon/pstpipeline/blob/main/docker#readme) for more details.

### I just want the data!

The raw data are rather large, so are shared here in the form of an ```R``` list saved as an ```.RDS``` file. See the [data-raw](https://github.com/qdercon/pstpipeline/tree/main/data-raw) folder and its accompanying [README](https://github.com/qdercon/pstpipeline/blob/main/data-raw#readme) for more details on how save the raw data and/or extract it as ```.csv``` files or ```pandas.DataFrame()``` objects.

## Key functions

- **parsing**:
    - ```import_single``` converts .txt JATOS outputs for a single participant to a list of cleaned dataframes.
    - ```import_multiple``` converts .txt JATOS outputs for multiple participants to a list of cleaned dataframes.


- **analysis and model checking**:
    - ```fit_learning_model``` automates using [CmdStanR](https://mc-stan.org/cmdstanr/) to run single and dual learning-rate $Q$-learning models [[4](#References)] in a hierarchical Bayesian manner, for both the training and test blocks, using either variational inference or MCMC. It is also capable of fitting a number of extended models which are additionally fit to trial-by-trial affect ratings [[2, 3](#References)]. Its functionality is comparable to [```hBayesDM::hBayesDM_model()```](https://github.com/CCS-Lab/hBayesDM/blob/develop/R/R/hBayesDM_model.R), but with a [CmdStan](https://github.com/stan-dev/cmdstan) backend.
    - ```generate_posterior_quantities``` enables posterior predictions for each MCMC sample to be generated in a separate session from a previously fitted model, via CmdStanR's [```$generate_quantities()```](https://mc-stan.org/cmdstanr/reference/model-method-generate-quantities.html) method.
    - ```parameter_glm``` is a wrapper for```cmdstan_glm```, a version of [```stan_glm```](https://mc-stan.org/rstanarm/reference/stan_glm.html) from the [```rstanarm```](https://github.com/stan-dev/rstanarm) package which has been modified to use [```cmdstanr```](https://mc-stan.org/cmdstanr/) rather than [```rstan```](https://cran.r-project.org/web/packages/rstan/index.html) as its backend (which greatly reduces the number of dependencies required). It is used to run adjusted Bayesian GLMs relating the individual-level posterior means of the learning parameters to outcome variable(s) of interest.
    - ```simulate_QL``` can be used to simulate data from the various $Q$-learing models, either using random samples from chosen distributions or the observed individual-level parameter estimates. The output can then be fit to the models, to check whether the parameters can be adequately recovered.


- **plotting**:
    - ```plot_factors``` produces a variety of plots to visually check the predictive accuracy of the lasso regression model used to predict transdiagnostic psychiatric symptom dimensions from a subset of questionnaire questions.
    - ```plot_import``` is a flexible function which enables the visualisation of various aspects of the observed data, including from the training and test phases of the PST, and the affect questions asked throughout the task. These can be presented for a single individual, or aggregated across all individuals (after applying exclusion criteria), and can also be used to compare groups based on any binary covariate.
    - ```check_learning_models``` is a simple function to output plots of the group-level means for each of the free parameters, plus some visual model checks for MCMC chains via the [```bayesplot```](https://mc-stan.org/bayesplot/) package (traces and rank histograms for each of the chains).
    - ```plot_affect``` plots the predictions from affect models against raw observations, by adjective (i.e., happiness, confidence, or engagemtent), either at a grouped- or individual-level. It can also plot the posterior mean weights obtained from this model. 
    - ```plot_ppc``` is a flexible plotting function to compare posterior predictions for both training and test data to their observed values, across participants.
    - ```plot_glm``` plots the results of Bayesian GLMs with both a boxplot (depicting HDIs), and a posterior density (```geom_half_violin```; adapted from the [RainCloudPlots](https://github.com/RainCloudPlots/RainCloudPlots) repository [[5](#References)]).
    - ```plot_recovery``` produces correlation plots for the observed and recovered $Q$-learning parameters (after running of ```fit_learning_model``` on simulated data from ```simulate_QL```), as well as confusion matrices.
    - ```plot_raincloud``` produces plots of the posterior mean parameter values or transdiagnostic factors, optionally by group.


- **other helper functions**
    - ```get_affect_ppc``` obtains individual-level posterior predictions and observed data for plotting from a affect model fit, plus returns model fit metrics ($R^2$, MAE, RMSE) for each individual. 
    - ```get_preds_by_chain``` automates the loading of posterior predictions obtained from running ```generate_posterior_quantities``` for plotting, by importing and summing the binary choice predictions chain-by-chain, and collating them into far smaller summaries. It also includes an optional method which can help prevent memory overload when loading large numbers of predictions.
        - ```get_subsample``` is a function to obtain a smaller subsample of individuals from the larger dataset, which may be helpful for demonstration purposes.
    - ```make_par_df``` creates a ```tibble::tibble()``` of model parameters by participant ID from the results of learning model fits.
    - ```quantile_hdi``` is a function to get arbitrary quantiles of a posterior distribution, based on [```HDIofMCMC```](https://github.com/CCS-Lab/hBayesDM/blob/develop/R/R/HDIofMCMC.R) from the [hBayesDM](https://github.com/CCS-Lab/hBayesDM) package.

## References

1.   W-Y. Ahn, N. Haines, L. Zhang, Revealing Neurocomputational Mechanisms of Reinforcement Learning and Decision-Making With the hBayesDM Package. *Comput. Psychiatry.* **1**, 24-57 (2017).

2.   R. B. Rutledge, N. Skandali, P. Dayan, R. J. Dolan. A computational and neural model of momentary subjective well-being. *Proc. Natl. Acad. Sci. U.S.A.* **111**(33), 12252-12257 (2014).

3. 	D. C. Jangraw, H. Keren, H. Sun, R. L. Bedder, R. B. Rutledge, F. Pereira, A. G. Thomas, D. S. Pine, C. Zheng, D. M. Nielson, A. Stringaris. A highly replicable decline in mood during rest and simple tasks. *Nat Hum Behav.* ***7***, 596–610 (2023).

4.   M. J. Frank, A. A. Moustafa, H. M. Haughey, T. Curran, K. E. Hutchison, Genetic triple dissociation reveals multiple roles for dopamine in reinforcement learning. *Proc. Natl. Acad. Sci. U.S.A.* **104**(41), 16311–16316 (2007).

5.   M. Allen, D. Poggiali, K. Whitaker et al. Raincloud plots: a multi-platform tool for robust data visualization [version 2; peer review: 2 approved]. *Wellcome Open Res.* **4**, 63 (2021).
