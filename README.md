# pstpipeline

[![R-CMD-check](https://github.com/qdercon/pstpipeline/actions/workflows/main.yml/badge.svg)](https://github.com/qdercon/pstpipeline/actions/workflows/main.yml)
[![Docker](https://img.shields.io/badge/dockerhub-image-important.svg?logo=Docker)](https://hub.docker.com/repository/docker/qdercon/pstpipeline/general#)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/qdercon/pstpipeline/blob/main/LICENSE)

**p**robabilistic **s**election **t**ask **pipeline**

### an R package to clean, analyse, and present data from a large online learning study

**Full methods and results from this study can be found in the preprint:**

> Q. Dercon<sup>†</sup>, S. Z. Mehrhof<sup>†</sup>, T. R. Sandhu, C. Hitchcock, R. P. Lawson, D. A. Pizzagalli, T. Dalgleish, C. L. Nord. A Core Component of Psychological Therapy Causes Adaptive Changes in Computational Learning Mechanisms. *PsyArXiv* (2022). https://psyarxiv.com/jmnek.

This package, which accompanies the preprint, is not meant to be a brand new toolkit - better, more flexible packages are available for analysing computational psychiatry experiments. Indeed, some of the R and most of the Stan code is inspired by and/or modified from other R packages, in particular [```hBayesDM```](https://github.com/CCS-Lab/hBayesDM) [[1](#References)] and [```rstanarm```](https://mc-stan.org/rstanarm/).

Instead, its main aims are as follows:

1.  To make it easier for our specific analyses to be replicated by others.
2.  To demonstrate a complete pre- and post-processing pipeline for a common learning task, which (hopefully) shows that such workflows are a) not overwhelmingly difficult to adopt, and b) can elicit valuable mechanistic insights.
3.  To do the above in a high-level manner, while still giving the user control over key aspects - most functionality of the package can be achieved with single-line function calls.

### Using the package

#### Docker container

The easiest way to interactively run all the analyses is to download and mount the Docker image. To do so, first download and install the relevant version of [Docker](https://docs.docker.com/get-docker/) for your OS, and then run the following in a command prompt:

```
docker pull qdercon/pstpipeline:v0.1.0
```

The image includes everything required to run the Jupyter notebooks both locally or on a cluster/cloud server in a containerised environment (i.e., local package installs etc. will not be affected). More specifically, it is a Linux environment containing:

* R v4.1.2 plus all package dependencies (see "DESCRIPTION" file for full details)
* Python v3.9.5 plus all dependencies, including rpy2 for running R code in Jupyter notebooks
* JupyterLab
* CmdStan v2.28.1
* Jupyter notebooks and required raw data

The image can also be built locally from the included Dockerfile. For example, on Windows, to clone the repo, extract the relevant data & notebooks, and build the Docker image, you could run the following:

```
git clone https://github.com/qdercon/pstpipeline
cd pstpipeline
cp -a data-raw docker/data-raw
cp -a notebooks docker/notebooks
cd docker
docker build -t pstpipeline-docker .
```

Once downloaded or built, to mount the image, run the following in a command prompt:

```
docker run -it --rm -p 8888:8888 -v <:/path/to/folder>:/root/<mount_folder_name>/ pstpipeline-docker
```

The -v flag and the path that follows is optional; this allows you to "mount" a folder on the disk to enable notebooks and model outputs to be saved locally. The command will output a link beginning with ```http//:127.0.0.1:8888/lab?token=``` which can be copied and pasted into a browser to open JupyterLab.

#### Local R installation

To install the R package and all dependencies directly, run the following:

```R
# install.packages("remotes")
remotes::install_github("qdercon/pstpipeline")
```

The majority of the code in the notebooks is in R, so if you wish to run things this way I would recommend taking a look at the notebooks, and then copy/paste or write your own function calls as appropriate. All user functions (listed below) are fully documented; this documentation can be easily accessed via the ```?``` prefix in R/RStudio. Though written primarily for our specific data/analyses, the functions are written to be relatively flexible, and aspects can be easily modified (e.g., to add new models).


### Key functions

- **parsing**:
    - ```import_single``` converts .txt JATOS outputs for a single participant to a list of cleaned dataframes.
    - ```import_multiple``` converts .txt JATOS outputs for multiple participants to a list of cleaned dataframes.


- **analysis and model checking**:
    - ```fit_learning_model``` automates using [CmdStanR](https://mc-stan.org/cmdstanr/) to run 1-alpha and 2-alpha Q-learning models [[2](#References)] in a hierarchical Bayesian manner, for both the training and test blocks, using either variational inference or MCMC. Its functionality is comparable to [```hBayesDM::hBayesDM_model()```](https://rdrr.io/cran/hBayesDM/src/R/hBayesDM_model.R), but with a [CmdStan](https://github.com/stan-dev/cmdstan) backend.
    - ```generate_posterior_quantities``` enables posterior predictions for each MCMC sample to be generated in a separate session from a previously fitted model, via CmdStanR's [```$generate_quantities()```](https://mc-stan.org/cmdstanr/reference/model-method-generate-quantities.html) method.
    - ```parameter_glm``` is a wrapper for```cmdstan_glm```, a version of [```stan_glm```](https://mc-stan.org/rstanarm/reference/stan_glm.html) from the [```rstanarm```](https://github.com/stan-dev/rstanarm) package which has been modified to use [```cmdstanr```](https://mc-stan.org/cmdstanr/) rather than [```rstan```](https://cran.r-project.org/web/packages/rstan/index.html) as its backend (which greatly reduces the number of dependencies required). It is used to run adjusted Bayesian GLMs relating the individual-level posterior means of the learning parameters to outcome variable(s) of interest.
    - ```simulate_QL``` can be used to simulate data from the various QL models, either using random samples from chosen distributions or the observed individual-level parameter estimates. The output can then be fit to the models, to check whether the parameters can be adequately recovered.


- **plotting**:
    - ```plot_factors``` produces a variety of plots to visually check the predictive accuracy of the lasso regression model used to predict transdiagnostic psychiatric symptom dimensions from a subset of questionnaire questions.
    - ```plot_import``` is a flexible function which enables the visualisation of various aspects of the observed data, including from the training and test phases of the PST, and the affect questions asked throughout the task. These can be presented for a single individual, or aggregated across all individuals (after applying exclusion criteria), and can also be used to compare groups based on any binary covariate.
    - ```check_learning_models``` is a simple function to output plots of the group-level means for each of the free parameters, plus some visual model checks for MCMC chains via the [```bayesplot```](https://mc-stan.org/bayesplot/) package (traces and rank histograms for each of the chains).
    - ```plot_ppc``` is a flexible plotting function to compare posterior predictions for both training and test data to their observed values, across participants.
    - ```plot_glm``` plots the results of Bayesian GLMs with both a boxplot (depicting HDIs), and a posterior density (```geom_half_violin```; adapted from the [RainCloudPlots](https://github.com/RainCloudPlots/RainCloudPlots) repository [[3](#References)]).
    - ```plot_recovery``` produces correlation plots for the observed and recovered QL parameters (after running of ```fit_learning_model``` on simulated data from ```simulate_QL```), as well as confusion matrices.
    - ```plot_raincloud``` produces plots of the posterior mean parameter values or transdiagnostic factors, optionally by group.


- **other helper functions**
    - ```get_subsample``` is a function to obtain a smaller subsample of individuals from the larger dataset, which may be helpful for demonstration purposes.
    - ```get_preds_by_chain``` automates the loading of posterior predictions obtained from running ```generate_posterior_quantities``` for plotting, by importing and summing the binary choice predictions chain-by-chain, and collating them into far smaller summaries. It also includes an optional method which can help prevent memory overload when loading large numbers of predictions.

### References

1.   W-Y. Ahn, N. Haines, L. Zhang, Revealing Neurocomputational Mechanisms of Reinforcement Learning and Decision-Making With the hBayesDM Package. *Comput. Psychiatry.* **1**, 24 (2017).

2.   M. J. Frank, A. A. Moustafa, H. M. Haughey, T. Curran, K. E. Hutchison, Genetic triple dissociation reveals multiple roles for dopamine in reinforcement learning. *Proc. Natl. Acad. Sci. U.S.A.* **104**, 16311–16316 (2007).

3.   M. Allen, D. Poggiali, K. Whitaker et al. Raincloud plots: a multi-platform tool for robust data visualization [version 2; peer review: 2 approved]. *Wellcome Open Res.* **4**, 63 (2021).
