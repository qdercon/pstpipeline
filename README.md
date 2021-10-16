# pstpipeline

**p**robabalistic **s**election **t**ask **pipeline**

#### an R package to clean, analyse, and present data from a large online learning study.

To install the R package and all dependencies:

```R
if (!require(remotes)) install.packages("remotes")
remotes::install_github("qdercon/pstpipeline")
```

Key elements of the package:

-   **parsing functions**:
    - ```import_single``` converts .txt JATOS outputs for a single participant to an R list 
    - ```import_multiple``` converts .txt JATOS outputs for multiple participants to an R list


-   **analysis functions**:
    - ```fit_learning_models``` automates using cmdstanr to run 1-alpha and 2-alpha Q-learning models [[1](#References)] in a hierarchical Bayesian manner, for both the training and test blocks.
    - ```generate_posterior_quantities``` enables posterior predictions and log-likelihoods for each MCMC sample to be generated in a separate session (as this is otherwise memory-intensive)


-   **model checks**:

-   **plotting functions**:

Please note much of the R and Stan code relies heavily on (and in some cases was adapted directly from) other packages, in particular [```hBayesDM```](https://github.com/CCS-Lab/hBayesDM) [[2](#References)]. This package is not meant to supercede or replace this work, nor will it be extensively updated to accomodate other workflows; its primary aim is make it easier for our specific analyses to be replicated by others.


#### References

1.   M. J. Frank, A. A. Moustafa, H. M. Haughey, T. Curran, K. E. Hutchison, Genetic triple dissociation reveals multiple roles for dopamine in reinforcement learning. *Proc. Natl. Acad. Sci. U.S.A.* **104**, 16311–16316 (2007).

2.   W-Y. Ahn, N. Haines, L. Zhang, Revealing Neurocomputational Mechanisms of Reinforcement Learning and Decision-Making With the hBayesDM Package. *Comput. Psychiatry.* **1**, 24 (2017).