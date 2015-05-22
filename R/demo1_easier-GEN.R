#' A cool function
#'
#' Well, not really cool. Just add 1 to x.
#' @param x a numeric vector
#' @export
#' @examples
#' add_one(1)
#' add_one(1:10)
add_one = function(x) {
  x + 1
}
#' MLE for the Gamma distribution
#'
#' Estimate the parameters (alpha and beta) of the Gamma distribution using
#' maximum likelihood.
#' @param data the data vector assumed to be generated from the Gamma
#'   distribution
#' @param start the initial values for the parameters of the Gamma distribution
#'   (passed to \code{\link{optim}()})
#' @param vcov whether to return an approximate variance-covariance matrix of
#'   the parameter vector
#' @return A list with elements \code{estimate} (parameter estimates for alpha
#'   and beta) and, if \code{vcov = TRUE}, \code{vcov} (the variance-covariance
#'   matrix of the parameter vector).
#' @export
mle_gamma = function(data, start = c(1, 1), vcov = FALSE) {
  loglike = function(param, x) {
    a = param[1]  # alpha (the shape parameter)
    b = param[2]  # beta (the rate parameter)
    n = length(x)
    n * (a * log(b) - lgamma(a)) + (a - 1) * sum(log(x)) - b * sum(x)
  }
  opt = optim(start, loglike, x = data, hessian = vcov, control = list(fnscale = -1))
  if (opt$convergence != 0) stop('optim() failed to converge')
  res = list(estimate = opt$par)
  if (vcov) res$vcov = solve(-opt$hessian)
  res
}
