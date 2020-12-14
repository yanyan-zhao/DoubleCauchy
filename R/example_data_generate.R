#' A function to simulate data for an example
#'
#' This function returns a list contains the simulated outcome variable Y and genotype matrix G.
#'
#' @param n sample size.
#' @param J dimension for genotype matrix.
#' @param beta a vector of parameters.
#'
#' @return A list contains outcome vector Y and genotype matrix G.
#'
#' @importFrom MASS mvrnorm
#' @import stats
#'
#' @examples
#' \dontrun{
#'  library(MASS)
#'  J = 100
#'  n = 50
#'  beta=c(rep(2,3),rep(0,J-3))
#'  data=example_data_generate(n, J, beta)
#' }
#' @export
#'

example_data_generate<-function(n, J, beta){

  #generate a AR(1) covariance matrix
  times <- 1:J
  del <- 0.5
  sigma <- 1
  H <- abs(outer(times, times, "-"))
  sigx <- sigma * del^H
  sigx[cbind(1:J, 1:J)] <- sigx[cbind(1:J, 1:J)] * sigma

  mu0=rep(0,J)
  G=mvrnorm(n,mu0,sigx)
  Gb = G%*%beta
  pr = exp(Gb) / ( 1 + exp(Gb))   ## inverse-logit
  Y = rbinom(n, 1, pr)

  results = list(Y = Y, G = G)
  return(results)
}






