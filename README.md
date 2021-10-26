# pstpipeline

## **p**robabilistic **s**election **t**ask **pipeline**

### an R package to clean, analyse, and present data from a large online learning study.

The idea of this package is not to be a brand new toolkit - better, more flexible packages are available for analysing computational psychiatry experiments. Indeed, much of the R and Stan code relies heavily on (and in some cases was adapted directly from) other packages, in particular [```hBayesDM```](https://github.com/CCS-Lab/hBayesDM) [[1](#References)]. This package is not meant to supercede or replace this work, nor will it be extensively updated to accomodate other workflows.

The aims of the package (and the interactive notebooks) are two-fold:

1.  To make it easier for our specific analyses to be replicated by others.
2.  To demonstrate a complete pre- and post-processing pipeline for a common learning task, which (hopefully) shows that such workflows are a) not overwhelmingly difficult to adopt, and b) can elicit valuable mechanistic insights.

#### 1. quickstart

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/qdercon/pstpipeline/)

The easiest way to interactively play around with the dataset and replicate the results is to open the Jupyter notebooks in Google Colabatory (i.e., click the button above). The majority of the code can be run completely on the cloud with no downloads necessary. 
The ```pstpipeline``` package and dependencies not already in Google Colab (most are) can be installed manually as explained in the notebooks. For speed, I would recommend sticking to variational inference if you are interested in (re)fitting models - it will elicit very similar results to MCMC and is orders of magnitude quicker. It also requires fewer R packages to be installed.

#### 2. running locally

To install the R package and all dependencies:

```R
if (!require(remotes)) install.packages("remotes")
remotes::install_github("qdercon/pstpipeline")
```

To get started, I would recommend taking a look at the Jupyter notebooks (under 'notebooks'). The vast majority of the code is R, so while these notebooks can easily be run locally using Jupyter Lab, this will additionally require a Python environment and the ryp2 library, and so it's probably easier to copy paste the relevant cells and run it in RStudio instead. Jupyter notebooks are included in the repo due to their ability to be opened in Google Colabatory, and because RMarkdown notebooks don't display nicely on GitHub.

#### 3. package components and use-cases

Key elements of the package:

- **parsing functions**:
    - ```import_single``` converts .txt JATOS outputs for a single participant to an R list 
    - ```import_multiple``` converts .txt JATOS outputs for multiple participants to an R list


- **analysis functions**:
    - ```fit_learning_models``` automates using [CmdStanR](https://mc-stan.org/cmdstanr/) to run 1-alpha and 2-alpha Q-learning models [[2](#References)] in a hierarchical Bayesian manner, for both the training and test blocks, using either variational inference or MCMC
    - ```generate_posterior_quantities``` enables posterior predictions and log-likelihoods for each MCMC sample to be generated in a separate session (as this is otherwise memory-intensive)


- **model checks**:
    - ```check_learning_models``` is a simple function to output plots of the group-level means for each of the free parameters, plus some visual model checks for MCMC chains (traces, rank histograms)


- **plotting functions**:
    - ```plot_factors``` is produces a variety of plots to visually check the predictive accuracy of the lasso regression model used to predict transdiagnostic psychiatric symptom dimensions from questionnaire questions
    - ```plot_import``` is a flexible function which enables the visualisation of various aspects of the observed data, including from the training and test phases of the PST, and the affect questions asked throughout the task. These can be presented for a single individual, or aggregated across all individuals (after applying exclusion criteria). Can also be used to compare groups based on any binary covariate    
    - ```plot_ppc``` is a flexible plotting function to compare posterior predictions for both training and test data to their observed values, across participants - can be plotted at an individual or group level


- **other helper functions**
    - ```get_subsample``` function to obtain a smaller subsample of individuals
    - ```save_preds_by_chain``` automates the memory-efficient loading of posterior predictions, likely obtained from running ```generate_posterior_quantities```, for plotting, by importing and summing the binary choice predictions chain-by-chain



### References

1.   W-Y. Ahn, N. Haines, L. Zhang, Revealing Neurocomputational Mechanisms of Reinforcement Learning and Decision-Making With the hBayesDM Package. *Comput. Psychiatry.* **1**, 24 (2017).

2.   M. J. Frank, A. A. Moustafa, H. M. Haughey, T. Curran, K. E. Hutchison, Genetic triple dissociation reveals multiple roles for dopamine in reinforcement learning. *Proc. Natl. Acad. Sci. U.S.A.* **104**, 16311–16316 (2007).
