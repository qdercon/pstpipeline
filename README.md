# pstpipeline

**p**robabilistic **s**election **t**ask **pipeline**

#### quickstart

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/qdercon/pstpipeline/)

#### an R package to clean, analyse, and present data from a large online learning study.

To install the R package and all dependencies:

```R
if (!require(remotes)) install.packages("remotes")
remotes::install_github("qdercon/pstpipeline")
```

Key elements of the package:

- **parsing functions**:
    - ```import_single``` converts .txt JATOS outputs for a single participant to an R list 
    - ```import_multiple``` converts .txt JATOS outputs for multiple participants to an R list


- **analysis functions**:
    - ```fit_learning_models``` automates using [CmdStanR](https://mc-stan.org/cmdstanr/) to run 1-alpha and 2-alpha Q-learning models [[1](#References)] in a hierarchical Bayesian manner, for both the training and test blocks, using either variational inference or MCMC
    - ```generate_posterior_quantities``` enables posterior predictions and log-likelihoods for each MCMC sample to be generated in a separate session (as this is otherwise memory-intensive)


- **model checks**:
    - ```check_learning_models``` is a simple function to output plots of the group-level means for each of the free parameters, plus some visual model checks for MCMC chains (traces, rank histograms)


- **plotting functions**:
    - ```plot_factors``` is produces a variety of plots to visually check the predictive accuracy of the lasso regression model used to predict transdiagnostic psychiatric symptom dimensions from questionnaire questions
    - ```plot_import``` is a flexible function which enables the visualisation of various aspects of the observed data, including from the training and test phases of the PST, and the affect questions asked throughout the task. These can be presented for a single individual, or aggregated across all individuals (after applying exclusion criteria). Can also be used to compare groups based on any binary covariate    
    - ```plot_PPCs``` is a flexible plotting function to compare posterior predictions for both training and test data to their observed values, across participants - can be plotted at an individual or group level


- **other helper functions**
    - ```get_subsample``` function to obtain a smaller subsample of individuals
    - ```load_PPCs``` automates the memory-efficient loading of posterior predictions, likely obtained from running ```generate_posterior_quantities```, by calculating means and SEs chain-by-chain. This is useful as, depending on how many MCMC samples selected, loading all predictions for all individuals at once can cause crashes


Please note much of the R and Stan code relies heavily on (and in some cases was adapted directly from) other packages, in particular [```hBayesDM```](https://github.com/CCS-Lab/hBayesDM) [[2](#References)]. This package is not meant to supercede or replace this work, nor will it be extensively updated to accomodate other workflows; its primary aim is make it easier for our specific analyses to be replicated by others.


#### References

1.   M. J. Frank, A. A. Moustafa, H. M. Haughey, T. Curran, K. E. Hutchison, Genetic triple dissociation reveals multiple roles for dopamine in reinforcement learning. *Proc. Natl. Acad. Sci. U.S.A.* **104**, 16311–16316 (2007).

2.   W-Y. Ahn, N. Haines, L. Zhang, Revealing Neurocomputational Mechanisms of Reinforcement Learning and Decision-Making With the hBayesDM Package. *Comput. Psychiatry.* **1**, 24 (2017).
