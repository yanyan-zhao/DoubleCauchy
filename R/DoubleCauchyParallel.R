#' A stable and adaptive polygenic signal detection based on repeated sample splitting.
#'
#' This function is a parallel version of DoubleCauchy function,
#' which can parallel the sample splitting times to improve
#' computational efficiency.
#'
#' @param n1 number of individuals for training.
#' @param m number of sample splitting.
#' @param Y The outcome variable of interest; a numeric vector of length n.
#' @param G a numeric n by J genotype matrix.
#' @param varselec.method variable selection methods for training sample.
#' @param J2 number of variables selected based on training sample. Only needed when \code{DCSIS} and \code{SIS} is selected.
#' @param alpha The elasticnet mixing parameter, with \eqn{0 \le \alpha \le 1}.
#'              \code{alpha=1} is the lasso penalty, and \code{alpha=0} the ridge penalty. Only needed when \code{ElasticNet} is selected.
#' @param cor.est The method for estimating a positive definite correlation matrix.
#' @param lam The tuning parameter for estimating correlation matrix; Could be either a scalar or
#'            a p by p symmetric matrix with an irrelevant diagonal while \code{pdsoft} is selected;
#'            Should be a vector when \code{pdsoft.cv} is selected.
#'            See \code{pdsoft} and \code{pdsoft.cv} in \code{PDSCE} for details.
#' @param pow.param an integer numeric vector specifying the \eqn{\gamma} values used in the test statistics.
#' @param seed seeds for sample splitting, should be a scalar.
#' @param family Specify response type. Only needed when \code{SIS} is selected and only "gaussian" and "binomial" are available.
#' @param ncores Number of cores to be used in parallel computing (default=1).
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
#' @importFrom SIS SIS
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
#'
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
#'  DoubleCauchyParallel(n1=25, m=10, Y=Y, G=G, varselec.method="SIS", J2=25,
#'  cor.est="pdsoft", lam=lam, pow.param=pow.param, family="binomial", ncores=2)
#'
#'
#' @keywords Parallel, Stable and adaptive test, Sample splitting, Double Cauchy
#'
#' @references Zhao and Sun (2020). A stable and adaptive polygenic signal detection method based on repeated sample splitting.
#'    \emph{arXiv:2008.02442}.
#'
#' @export
#'
DoubleCauchyParallel<-function(n1, m, Y, G, varselec.method=c("DCSIS","ElasticNet","SIS"),
                       J2=NULL, alpha=NULL, cor.est=c("pdsoft","pdsoft.cv"),
                       lam, pow.param=c(0:10), seed=NULL,family=c("gaussian","binomial"), ncores=1){

  n=dim(G)[1]
  J=dim(G)[2]

  index=c(rep(1,n1),rep(0,n-n1))

  varmeth <- match.arg(varselec.method)
  cor.est <- match.arg(cor.est)
  family <- match.arg(family)

  pfinal=NULL
  i=NULL
  if(is.null(seed)){
    ##################################
    #  parallel
    ##################################
    registerDoParallel(ncores)
    pfinal<-foreach(i=1:m, .combine=rbind, .errorhandling='pass')%dopar%{
      set.seed(i)
      ind_random=as.logical(sample(index))
      samsplit(Y=Y, G=G, ind_random=ind_random, varselec.method=varselec.method, J2=J2, alpha=alpha, cor.est=cor.est, lam=lam, pow.param=pow.param,family=family)
    }
   stopImplicitCluster()

  }else{
    ##################################
    #  parallel
    ##################################
    registerDoParallel(ncores)
    pfinal<-foreach(i=1:m, .combine=rbind, .errorhandling='pass')%dopar%{
      set.seed((seed-1)*m+i)
      ind_random=as.logical(sample(index))
      samsplit(Y=Y, G=G, ind_random=ind_random, varselec.method=varselec.method, J2=J2, alpha=alpha, cor.est=cor.est, lam=lam, pow.param=pow.param,family=family)
    }
    stopImplicitCluster()

  }

  cauchyfun<-function(x){
    t0=(1/length(x))*sum(tan((0.5-x)*pi))
    p_cauchy=0.5-(atan(t0)/pi)
    return(p_cauchy)
  }

  p.results=apply(pfinal, 2, cauchyfun)
  return(p.results)
}







