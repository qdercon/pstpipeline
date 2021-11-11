# Adapted from the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016, 2017 Trustees of Columbia University
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#' Bayesian generalized linear models via CmdStan
#'
#' Generalized linear modeling for Gaussian and gamma responses, with optional prior
#' distributions for the coefficients, intercept, and auxiliary parameters.
#'
#' @param formula,data,subset Same as \code{\link[stats]{glm}},
#' but \emph{we strongly advise against omitting the \code{data}
#' argument}.
#' @param family Same as \code{\link[stats]{glm}}. Accepted families for this function
#' are \code{gaussian()} or \code{Gamma()}, and any link function.
#' @param model,offset,weights Same as \code{\link[stats]{glm}}.
#' @param na.action,contrasts Same as \code{\link[stats]{glm}}, but
#' rarely specified.
#' @param y Logical scalar indicating whether to return the response vector.
#' @param x Logical scalar indicating whether to return the design matrix.
#' @param algorithm Argument "sampling" is for MCMC (default), while "meanfield"
#' and "fullrank" are variational algorithms ("meanfield" is the CmdStan default).
#' @param out_dir Output directory for model fit environment.
#' @param ... Further arguments passed to [cmdstanr::sample] (i.e., \code{refresh},
#'   \code{iter_warmup}, \code{iter_sampling}, \code{chains} etc.)
#' @param prior The prior distribution for the (non-hierarchical) regression
#'   coefficients.
#'
#'   The default priors are described in the vignette
#'   \href{http://mc-stan.org/rstanarm/articles/priors.html}{\emph{Prior
#'   Distributions for rstanarm Models}}.
#'   If not using the default, \code{prior} should be a call to one of the
#'   various functions provided by \pkg{rstanarm} for specifying priors. The
#'   subset of these functions that can be used for the prior on the
#'   coefficients can be grouped into several "families":
#'
#'   \tabular{ll}{
#'     \strong{Family} \tab \strong{Functions} \cr
#'     \emph{Student t family} \tab \code{normal}, \code{student_t}, \code{cauchy} \cr
#'     \emph{Hierarchical shrinkage family} \tab \code{hs}, \code{hs_plus} \cr
#'     \emph{Laplace family} \tab \code{laplace}, \code{lasso} \cr
#'     \emph{Product normal family} \tab \code{product_normal} \cr
#'   }
#'
#'   See \[http://mc-stan.org/rstanarm/reference/priors.html](here) for details on the
#'   families and how to specify the arguments for all of the functions in the table above.
#'   To omit a prior ---i.e., to use a flat (improper) uniform prior---
#'   \code{prior} can be set to \code{NULL}, although this is rarely a good
#'   idea.
#'
#'   \strong{Note:} Unless \code{QR=TRUE}, if \code{prior} is from the Student t
#'   family or Laplace family, and if the \code{autoscale} argument to the
#'   function used to specify the prior is left at its default and recommended
#'   value of \code{TRUE}, then the default or user-specified prior scale(s) may
#'   be adjusted internally based on the scales of the predictors.
#'
#' @param prior_intercept The prior distribution for the intercept (after
#'   centering all predictors, see note below).
#'
#'   The default prior is described in the vignette
#'   \href{http://mc-stan.org/rstanarm/articles/priors.html}{\emph{Prior
#'   Distributions for rstanarm Models}}.
#'   If not using the default, \code{prior_intercept} can be a call to
#'   \code{normal}, \code{student_t} or \code{cauchy}. To omit a
#'   prior on the intercept ---i.e., to use a flat (improper) uniform prior---
#'   \code{prior_intercept} can be set to \code{NULL}.
#'
#'   \strong{Note:} If using a dense representation of the design matrix
#'   ---i.e., if the \code{sparse} argument is left at its default value of
#'   \code{FALSE}--- then the prior distribution for the intercept is set so it
#'   applies to the value \emph{when all predictors are centered} (you don't
#'   need to manually center them). This is explained further in
#'   [Prior Distributions for rstanarm Models](https://mc-stan.org/rstanarm/articles/priors.html)
#'   If you prefer to specify a prior on the intercept without the predictors
#'   being auto-centered, then you have to omit the intercept from the
#'   \code{\link[stats]{formula}} and include a column of ones as a predictor,
#'   in which case some element of \code{prior} specifies the prior on it,
#'   rather than \code{prior_intercept}. Regardless of how
#'   \code{prior_intercept} is specified, the reported \emph{estimates} of the
#'   intercept always correspond to a parameterization without centered
#'   predictors (i.e., same as in \code{glm}).
#'
#' @param prior_aux The prior distribution for the "auxiliary" parameter (if
#'   applicable). The "auxiliary" parameter refers to a different parameter
#'   depending on the \code{family}. For Gaussian models \code{prior_aux}
#'   controls \code{"sigma"}, the error
#'   standard deviation. For negative binomial models \code{prior_aux} controls
#'   \code{"reciprocal_dispersion"}, which is similar to the
#'   \code{"size"} parameter of \code{\link[stats:NegBinomial]{rnbinom}}:
#'   smaller values of \code{"reciprocal_dispersion"} correspond to
#'   greater dispersion. For gamma models \code{prior_aux} sets the prior on
#'   to the \code{"shape"} parameter (see e.g.,
#'   \code{\link[stats:GammaDist]{rgamma}}), and for inverse-Gaussian models it is the
#'   so-called \code{"lambda"} parameter (which is essentially the reciprocal of
#'   a scale parameter). Binomial and Poisson models do not have auxiliary
#'   parameters.
#'
#'   The default prior is described in the vignette
#'   \href{http://mc-stan.org/rstanarm/articles/priors.html}{\emph{Prior
#'   Distributions for rstanarm Models}}.
#'   If not using the default, \code{prior_aux} can be a call to
#'   \code{exponential} to use an exponential distribution, or \code{normal},
#'   \code{student_t} or \code{cauchy}, which results in a half-normal, half-t,
#'   or half-Cauchy prior. See [rstanarm::priors()] for details on these
#'   functions. To omit a prior ---i.e., to use a flat (improper) uniform
#'   prior--- set \code{prior_aux} to \code{NULL}.
#' @param prior_PD A logical scalar (defaulting to \code{FALSE}) indicating
#'   whether to draw from the prior predictive distribution instead of
#'   conditioning on the outcome.
#' @param mean_PPD A logical value indicating whether the sample mean of the
#'   posterior predictive distribution of the outcome should be calculated in
#'   the \code{generated quantities} block. If \code{TRUE} then \code{mean_PPD}
#'   is computed and displayed as a diagnostic in the printed output. A useful
#'   heuristic is to check if \code{mean_PPD} is plausible when compared to
#'   \code{mean(y)}. If it is plausible then this does \emph{not} mean that the
#'   model is good in general (only that it can reproduce the sample mean), but
#'   if \code{mean_PPD} is implausible then there may be something wrong, e.g.,
#'   severe model misspecification, problems with the data and/or priors,
#'   computational issues, etc.
#' @param sparse A logical scalar (defaulting to \code{FALSE}) indicating
#'   whether to use a sparse representation of the design (X) matrix.
#'   If \code{TRUE}, the the design matrix is not centered (since that would
#'   destroy the sparsity) and likewise it is not possible to specify both
#'   \code{QR = TRUE} and \code{sparse = TRUE}. Depending on how many zeros
#'   there are in the design matrix, setting \code{sparse = TRUE} may make
#'   the code run faster and can consume much less RAM.
#'
#' @details The \code{stan_glm} function is similar in syntax to
#'   \code{\link[stats]{glm}} but rather than performing maximum likelihood
#'   estimation of generalized linear models, full Bayesian estimation is
#'   performed (if \code{algorithm} is \code{"sampling"}) via MCMC. The Bayesian
#'   model adds priors (independent by default) on the coefficients of the GLM.
#'
#' @return A \code{cmdstanr::CmdStanMCMC()} object.
#' @export

cmdstan_glm <-
  function(formula,
           family = gaussian(),
           data,
           weights,
           subset,
           na.action = NULL,
           offset = NULL,
           model = TRUE,
           algorithm = c("sampling", "meanfield", "fullrank"),
           x = FALSE,
           y = TRUE,
           contrasts = NULL,
           out_dir = NULL,
           ...,
           prior = default_prior_coef(family),
           prior_intercept = default_prior_intercept(family),
           prior_aux = exponential(autoscale=TRUE),
           prior_PD = FALSE,
           mean_PPD = !prior_PD,
           sparse = FALSE) {

  family <- validate_family(family)
  data <- validate_data(data, if_missing = environment(formula))

  call <- match.call(expand.dots = TRUE)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "subset", "weights", "na.action", "offset"),
             table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$data <- data
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mf <- check_constant_vars(mf)
  mt <- attr(mf, "terms")
  Y <- array1D_check(model.response(mf, type = "any"))
  if (is.empty.model(mt))
    stop("No intercept or predictors specified.", call. = FALSE)
  X <- model.matrix(mt, mf, contrasts)
  contrasts <- attr(X, "contrasts")
  weights <- validate_weights(as.vector(model.weights(mf)))
  offset <- validate_offset(as.vector(model.offset(mf)), y = Y)

  if (prior_PD) {
    # can result in errors if draws from prior are weird
    mean_PPD <- FALSE
  }

  stanfit <- cmdstan_glm.fit(
    x = X,
    y = Y,
    weights = weights,
    offset = offset,
    family = family,
    algorithm = algorithm,
    out_dir = out_dir,
    prior = prior,
    prior_intercept = prior_intercept,
    prior_aux = prior_aux,
    prior_PD = prior_PD,
    mean_PPD = mean_PPD,
    sparse = sparse,
    ...
  )

  return(stanfit)
}


# rstanarm::stan_glm internal fns ---------------------------------------------------------

validate_weights <- function(w) {
  if (missing(w) || is.null(w)) {
    w <- double(0)
  } else {
    if (!is.numeric(w))
      stop("'weights' must be a numeric vector.",
           call. = FALSE)
    if (any(w < 0))
      stop("Negative weights are not allowed.",
           call. = FALSE)
  }

  return(w)
}

validate_offset <- function(o, y) {
  if (is.null(o)) {
    o <- double(0)
  } else {
    if (length(o) != NROW(y))
      stop(gettextf("Number of offsets is %d but should be %d (number of observations)",
                    length(o), NROW(y)), domain = NA, call. = FALSE)
  }
  return(o)
}
