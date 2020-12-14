#' An adaptive test leveraging side information
#'
#' The function returns the p values of an adaptive test by leveraging available side information.
#'
#' @param Y The outcome variable of interest; a numeric vector of length n.
#' @param G a numeric n by J genotype matrix.
#' @param cor.est The method for estimating a positive definite correlation matrix.
#' @param lam The tuning parameter for estimating correlation matrix; Could be either a scalar or
#'            a p by p symmetric matrix with an irrelevant diagonal while \code{pdsoft} is selected;
#'            Should be a vector when \code{pdsoft.cv} is selected.
#'            See \code{pdsoft} and \code{pdsoft.cv} in \code{PDSCE} for details.
#' @param weights a numeric vector of length J, contains the weights from side information.
#' @param pow.param an integer numeric vector specifying the \eqn{\gamma} values used in the test statistics.
#'
#' @return A vector of p values.
#'
#' @importFrom matrixcalc is.positive.definite
#' @importFrom MASS mvrnorm
#' @importFrom CompQuadForm davies
#' @importFrom PDSCE pdsoft pdsoft.cv
#' @importFrom VariableScreening screenIID
#' @import stats
#'
#' @examples
#'  J = 100
#'  n = 50
#'  beta=c(rep(2,3),rep(0,J-3))
#'  data=example_data_generate(n, J, beta)
#'  Y=data$Y
#'  G=data$G
#'  lam=0.1
#'  weights=c(J:1)
#'  pow.param=c(1:5)
#'  AdapSide(Y, G, cor.est="pdsoft", lam, weights, pow.param)
#'
#'
#' @keywords Side information, Adaptive test
#'
#' @references Zhao and Sun (2020). A stable and adaptive polygenic signal detection method based on repeated sample splitting.
#'    \emph{arXiv:2008.02442}.
#'
#' @export
#'

AdapSide=function(Y, G, cor.est=c("pdsoft","pdsoft.cv"), lam, weights, pow.param=c(0:10)){

  cor.est <- match.arg(cor.est)

  n=dim(G)[1]
  J=dim(G)[2]

  if (!is.vector(Y)) {
    stop("Y should be a numeric vector")
  }

  if (!is.matrix(G)) {
    stop("G should be a numeric matrix")
  }

  if (length(weights)!=J) {
    stop("Please provide a numerical vector of length equals to the dimension of G")
  }

  ######################
  # generate score vector
  ######################
  S=rep(0,J)

  for(j in 1:J){
    x=scale(G[,j])
    x=x[!is.na(x)]
    y=scale(Y)[!is.na(x)]
    nx=length(x)
    S[j]=as.numeric((t(y)%*%x)/sqrt(nx-1))
  }

  cov_G=cov(G)
  pd.test=is.positive.definite(cov_G, tol=1e-8)

  if(pd.test){
    sig0=cor(G)
  }else if(cor.est=="pdsoft"){
    output=pdsoft(cov_G, lam=lam, standard = T)
    sig0=output$theta
  }else if(cor.est=="pdsoft.cv"){
    output=pdsoft.cv(G, lam.vec=lam, standard = T, nsplits=5)
    corr_pd=output$sigma
    sig0=cov2cor(corr_pd)
  }else{
    stop("The estimation of correlation matrix of G is not positive definite, please specify a method for cor.est")
  }

  A=chol(sig0)
  p.vec=rep(NA,length(pow.param))

  for(i in 1:length(pow.param)){
      R=diag(weights^pow.param[i])
      D=A%*%R%*%t(A)
      egi.val=eigen(D)$values

      Tgamma= as.numeric(S%*%R%*%S)
      p.gamma=davies(Tgamma,lambda=egi.val)$Qq
      p.vec[i]=p.gamma
  }

  t0=(1/length(p.vec))*sum(tan((0.5-p.vec)*pi))
  p.cauchy=0.5-(atan(t0)/pi)

  pval=c(p.vec, p.cauchy)
  names(pval) <- c(paste("p.", pow.param, sep = ""), "p.cauchy")
  return(pval)
}







