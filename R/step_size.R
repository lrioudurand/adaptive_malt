
# setwd("~/Code/adaptive_malt/R")
source("pseudo_code.R")

#'Init stepszie
#' @description Heuristics for finding good initial step size [1, Algorithm 4]
#' [1]: Hoffman and Gelman, 2014, No-U Turn Sampler.
#' 
#' @param U A potential function to return the log-density of the distribution to be sampled from, up to an additive constant. It should input a real vector of the same length as \code{init} and output a scalar.
#' @param grad A function to return the gradient of the potential. It should input and output a real vector of the same length as \code{init}.
#' @param pos0 Initial position
#' @param tau0 Initial trajectory length
#' @param h0 Initial time step.
#' @param w Initial eigen vector.

init_stepsize = function(U, grad, pos0, tau0, w0, h0 = 1) {
  m_dummy = rep(0, length(pos0))
  Delta = osam(U, grad, pos0, tau0, h0, w0, m_dummy, iter = 1)$Delta
  a = 2 * (Delta < log(2)) - 1  # check exp(- Delta) > 0.5
  while (a * Delta < a * log(2)) {
    h0 = 2^a * h0
    Delta = osam(U, grad, pos0, tau0, h0, w0, m_dummy, iter = 1)$Delta
  }
  
  return(h0)
}
