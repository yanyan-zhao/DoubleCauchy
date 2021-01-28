samsplit=function(Y, G, ind_random, varselec.method=c("DCSIS","ElasticNet","SIS"),
                  J2=NULL, alpha=NULL, cor.est=c("pdsoft","pdsoft.cv"),
                  lam, pow.param=c(0:10),family=c("gaussian","binomial")){

  varmeth <- match.arg(varselec.method)
  cor.est <- match.arg(cor.est)
  pow.param <- pow.param
  family <- match.arg(family)

  if (!is.vector(Y)) {
    stop("Y should be a numeric vector")
  }

  if (!is.matrix(G)) {
    stop("G should be a numeric matrix")
  }

  y.train<-scale(Y[ind_random])
  G.train=scale(G[ind_random,])

  ind_random_2=as.logical(1-ind_random)
  y.test=scale(Y[ind_random_2])
  G.test=scale(G[ind_random_2,])

  if(varmeth == "SIS"){

    if (is.null(J2)) {
      stop("Please specify J2 to be used for variable selection in SIS")
    }

    Sresult = SIS(x = G.train, y = y.train, family=family, nsis=J2, iter=F, penalty="lasso")
    G.train.selec=G.train[,Sresult$sis.ix0]
    G.test.selec=G.test[,Sresult$sis.ix0]

    num.var=dim(G.test.selec)[2]

  }else if(varmeth == "DCSIS"){

    if (is.null(J2)) {
      stop("Please specify J2 to be used for variable selection in DCSIS")
    }

    DCres = screenIID(X = G.train, Y = y.train, method="DC-SIS")
    G.train.selec=G.train[,DCres$rank<=J2]
    G.test.selec=G.test[,DCres$rank<=J2]

    num.var=dim(G.test.selec)[2]

  }else{

    if (is.null(alpha)) {
      stop("Please specify alpha to be used for variable selection in ElasticNet")
    }

    cvfit=cv.glmnet(G.train, y.train, alpha=alpha)
    coef_s1=coef(cvfit, s="lambda.min")
    pp=1:dim(coef_s1)[1]
    variable_set=pp[coef_s1[,1]!= 0] ### returns nonzero coefs
    variable_set=variable_set[-1]
    variable_set=variable_set-1
    num.var=length(variable_set)

    if (num.var!=0){
      G.train.selec=G.train[,variable_set]
      G.test.selec=G.test[,variable_set]
    }else{
      stop("No variables are selected by ElasticNet, try specify a smaller alpha to relax the penalty")
    }
  }

  if(num.var!=0){
    # generate weights
    weights=rep(0,num.var)

    for(j in 1:num.var){
      x=G.train.selec[,j]
      nx=length(x)
      weights[j]=as.numeric((t(y.train)%*%x)/sqrt(nx-1))
    }

    y.test=as.vector(y.test)
    G.test.selec=as.matrix(G.test.selec)

    results=AdapSide(Y=y.test, G=G.test.selec, cor.est=cor.est, lam=lam, weights=weights, pow.param=pow.param)
    return(results)

  }else{
    stop("No variables are selected")
  }
}
