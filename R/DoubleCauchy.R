#' A stable and adaptive polygenic signal detection based on repeated sample splitting.
#'
#' The function returns the p values of a stable and adaptive test for
#' simultaneously testing regression coefficients
#' of generalized linear models of high dimensional data.
#'
#'
#'
#' @param n1 number of individuals for training.
#' @param m number of sample splitting.
#' @param Y The outcome variable of interest; a numeric vector of length n.
#' @param G a numeric n by J genotype matrix.
#' @param varselec.method variable selection methods for training sample.
#' @param J2 number of variables selected based on training sample. Only needed when \code{DCSIS} is selected.
#' @param alpha The elasticnet mixing parameter, with \eqn{0\le\alpha\le 1}.
#'              \code{alpha=1} is the lasso penalty, and \code{alpha=0} the ridge penalty. Only needed when \code{glmnet} is selected.
#' @param cor.est The method for estimating a positive definite correlation matrix.
#' @param lam The tuning parameter for estimating correlation matrix; Could be either a scalar or
#'            a p by p symmetric matrix with an irrelevant diagonal while \code{pdsoft} is selected;
#'            Should be a vector when \code{pdsoft.cv} is selected.
#'            See \code{pdsoft} and \code{pdsoft.cv} in \code{PDSCE} for details.
#' @param pow.param an integer numeric vector specifying the \eqn{\gamma} values used in the test statistics.
#' @param seed seeds for sample splitting, should be a scalar.
#'
#' @return A vector of p values.
#'
#' @importFrom matrixcalc is.positive.definite
#' @importFrom MASS mvrnorm
#' @importFrom CompQuadForm davies
#' @importFrom PDSCE pdsoft pdsoft.cv
#' @importFrom glmnet cv.glmnet
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
#'  n1=20
#'  m=10
#'  DoubleCauchy(n1, m, Y, G, varselec.method="DCSIS", J2=10,
#'  cor.est="pdsoft", lam=lam, pow.param=pow.param)
#'
#'
#' @keywords Stable and adaptive test, Sample splitting, Double Cauchy
#'
#' @references Zhao and Sun (2020). A stable and adaptive polygenic signal detection method based on repeated sample splitting.
#'    \emph{arXiv:2008.02442}.
#'
#' @export
#'
DoubleCauchy<-function(n1, m, Y, G, varselec.method=c("DCSIS","ElasticNet"),
                       J2=NULL, alpha=NULL, cor.est=c("pdsoft","pdsoft.cv"),
                       lam, pow.param=c(0:10), seed=NULL){

  n=dim(G)[1]
  J=dim(G)[2]

  index=c(rep(1,n1),rep(0,n-n1))

  varmeth <- match.arg(varselec.method)
  cor.est <- match.arg(cor.est)

  pfinal=NULL
  if(is.null(seed)){
    for(i in 1:m){
      set.seed(i)
      ind_random=as.logical(sample(index))
      pfinal0=samsplit(Y=Y, G=G, ind_random=ind_random, varselec.method=varselec.method, J2=J2, alpha=alpha, cor.est=cor.est, lam=lam, pow.param=pow.param)
      pfinal=rbind(pfinal,pfinal0)
    }
  }else{
    for(i in 1:m){
      set.seed(i*seed)
      ind_random=as.logical(sample(index))
      pfinal0=samsplit(Y=Y, G=G, ind_random=ind_random, varselec.method=varselec.method, J2=J2, alpha=alpha, cor.est=cor.est, lam=lam, pow.param=pow.param)
      pfinal=rbind(pfinal,pfinal0)
    }
  }

  cauchyfun<-function(x){
    t0=(1/length(x))*sum(tan((0.5-x)*pi))
    p_cauchy=0.5-(atan(t0)/pi)
    return(p_cauchy)
  }

  p.results=apply(pfinal, 2, cauchyfun)
  return(p.results)
}







