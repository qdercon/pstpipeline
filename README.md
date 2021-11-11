# pstpipeline

![R CMD check](https://github.com/qdercon/pstpipeline/actions/workflows/main.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://cran.r-project.org/web/licenses/MIT)

**p**robabilistic **s**election **t**ask **pipeline**

### an R package to clean, analyse, and present data from a large online learning study.

The idea of this package is not to be a brand new toolkit - better, more flexible packages are available for analysing computational psychiatry experiments. Indeed, much of the R and Stan code relies heavily on code modified from other packages, in particular [```hBayesDM```](https://github.com/CCS-Lab/hBayesDM) [[1](#References)], which inspired both the approach and execution of the analysis itself. 

This package is not meant to supercede or replace this work; instead its main aims are as follows:

1.  To make it easier for our specific analyses to be replicated by others.
2.  To demonstrate a complete pre- and post-processing pipeline for a common learning task, which (hopefully) shows that such workflows are a) not overwhelmingly difficult to adopt, and b) can elicit valuable mechanistic insights.
3.  To do the above in a high-level manner, while still giving the user control over key aspects - most functionality of the package can be achieved with single-line function calls.

#### Using the package

The easiest way to interactively run all the analyses is to download and mount the Docker image. This image contains everything required to run the Jupyter notebooks both locally or even on a cloud server (e.g., Google Cloud) in a containerised environment (i.e., local package installs etc. will not be affected).

Specifically, it contains:

* All R package dependencies (see "DESCRIPTION" file for full details)
* Python dependencies, including rpy2 for running R code in Jupyter notebooks
* Jupyter Lab
* CmdStan v2.28.1

To mount the image, and open a Jupyter notebook in your browser, run the following:

```
docker run -it --rm -p 8888:8888 -v [:/Path/To/Folder]:/root/local_mount/ pstpipeline-docker:latest
```

The -v flag and the path that follows is optional; this allows you to mount a local drive to enable notebooks/model outputs to be saved to locally should you wish.

Alternatively, to install the R package and all dependencies locally, run the following:

```R
# install.packages("remotes")
remotes::install_github("qdercon/pstpipeline")
```

The majority of the code in the notebooks is in R, so if you wish to run things this way I would recommend taking a look at the notebooks/R documentation, and then copy/paste or write your own function calls as appropriate.


#### Key components of the package

- **parsing functions**:
    - ```import_single``` converts .txt JATOS outputs for a single participant to a list of cleaned dataframes.
    - ```import_multiple``` converts .txt JATOS outputs for multiple participants to a list of cleaned dataframes.


- **analysis and model checking**:
    - ```fit_learning_model``` automates using [CmdStanR](https://mc-stan.org/cmdstanr/) to run 1-alpha and 2-alpha Q-learning models [[2](#References)] in a hierarchical Bayesian manner, for both the training and test blocks, using either variational inference or MCMC.
    - ```generate_posterior_quantities``` enables posterior predictions for each MCMC sample to be generated in a separate session using a previously fitted model (as this can be memory-intensive otherwise).
    - ```parameter_glm``` is a wrapper for```cmdstan_glm```, a version of [```stan_glm```](https://mc-stan.org/rstanarm/reference/stan_glm.html) from the [```rstanarm```](https://github.com/stan-dev/rstanarm) package which has been modified to use [```cmdstanr```](https://mc-stan.org/cmdstanr/) rather than [```rstan```](https://cran.r-project.org/web/packages/rstan/index.html) as its backend (which greatly reduces the number of dependencies required). It is used to run adjusted Bayesian GLMs relating the individual-level posterior means of the learning parameters to outcome variable(s) of interest.
    - ```simulate_QL``` can be used to simulate data from the various QL models, either using random samples from chosen distributions or the observed individual-level parameter estimates. The output can then be fit to the models, to check whether the parameters can be adequately recovered.


- **plotting functions**:
    - ```plot_factors``` produces a variety of plots to visually check the predictive accuracy of the lasso regression model used to predict transdiagnostic psychiatric symptom dimensions from a subset of questionnaire questions.
    - ```plot_import``` is a flexible function which enables the visualisation of various aspects of the observed data, including from the training and test phases of the PST, and the affect questions asked throughout the task. These can be presented for a single individual, or aggregated across all individuals (after applying exclusion criteria), and can also be used to compare groups based on any binary covariate.
    - ```check_learning_models``` is a simple function to output plots of the group-level means for each of the free parameters, plus some visual model checks for MCMC chains (traces and rank histograms for each of the chains).
    - ```plot_ppc``` is a flexible plotting function to compare posterior predictions for both training and test data to their observed values, across participants.
    - ```plot_glm``` plots the results of Bayesian GLMs with both a boxplot (depicting HDIs), and a posterior density made up of the posterior draws themselves, using a modified version of [```geom_quasirandom```](https://www.rdocumentation.org/packages/ggbeeswarm/versions/0.5.3/topics/geom_quasirandom) from the [```ggbeeswarm```](https://github.com/eclarke/ggbeeswarm) package.
    - ```plot_recovery``` produces correlation plots for the observed and recovered QL parameters (after running of ```fit_learning_model``` on simulated data from ```simulate_QL```), as well as confusion matrices.


- **other helper functions**
    - ```get_subsample``` is a function to obtain a smaller subsample of individuals from the larger dataset, which may be helpful for demonstration purposes.
    - ```get_preds_by_chain``` automates the loading of posterior predictions obtained from running ```generate_posterior_quantities``` for plotting, by importing and summing the binary choice predictions chain-by-chain, and collating them into far smaller summaries. It also includes an optional method which can help prevent memory overload when loading large numbers of predictions.
    


### References

1.   W-Y. Ahn, N. Haines, L. Zhang, Revealing Neurocomputational Mechanisms of Reinforcement Learning and Decision-Making With the hBayesDM Package. *Comput. Psychiatry.* **1**, 24 (2017).

2.   M. J. Frank, A. A. Moustafa, H. M. Haughey, T. Curran, K. E. Hutchison, Genetic triple dissociation reveals multiple roles for dopamine in reinforcement learning. *Proc. Natl. Acad. Sci. U.S.A.* **104**, 16311–16316 (2007).
