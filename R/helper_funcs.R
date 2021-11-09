#' @noRd
#' @keywords internal

single_hdi <- function(vars, cred) {
  # adapted from hBayesDM::HDIofMCMC

  sampleVec <- as.vector(t(vars))
  sortedPts = sort(sampleVec)
  ciIdxInc = floor(cred * length(sortedPts))
  nCIs = length(sortedPts) - ciIdxInc
  ciWidth = rep(0 , nCIs)
  for (i in 1:nCIs) {
    ciWidth[i] = sortedPts[i + ciIdxInc] - sortedPts[i]
  }
  HDImin = sortedPts[which.min(ciWidth)]
  HDImax = sortedPts[which.min(ciWidth) + ciIdxInc]
  HDIlim = c(HDImin , HDImax)
  return(as.vector(t(HDIlim)))
}
quantile_hdi <- function(var, quantile, transform = FALSE) {

  if (transform) {
   var <- log(var/100 + 1)
  }

  returns <- vector(mode = "numeric")
  for (q in quantile) {
    if (q == 0.5) {
      returns <- cbind(returns, median(var))
    } else if (q == 0) {
      returns <- cbind(returns, min(var))
    } else if (q == 1){
      returns <- cbind(returns, max(var))
    } else if (q < 0.5) {
      cred_mass = 1 - 2*q
      HDI_lower <- single_hdi(vars = var, cred = cred_mass)[1]
      returns <- cbind(returns, HDI_lower)
    } else {
      cred_mass = 2*q - 1
      HDI_upper <- single_hdi(vars = var, cred = cred_mass)[2]
      returns <- cbind(returns, HDI_upper)
    }
  }

  ret <- sort(returns)
  if (transform) {
    ret <- (exp(ret) - 1) *100
  }
  names(ret) <- sapply(1:length(quantile), FUN = function (x) paste0(quantile[x] * 100, "%"))

  return(ret)
}
family_ch <- function(param) {
  if (grepl("alpha", param)) return(Gamma(link = "log"))
  else return(gaussian())
}
