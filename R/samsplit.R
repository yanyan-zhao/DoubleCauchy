samsplit=function(Y, G, ind_random, varselec.method=c("DCSIS","ElasticNet"),
                  J2=NULL, alpha=NULL, cor.est=c("pdsoft","pdsoft.cv"),
                  lam, pow.param=c(0:10)){

  varmeth <- match.arg(varselec.method)
  cor.est <- match.arg(cor.est)
  pow.param <- pow.param

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

  if(varmeth == "DCSIS"){

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
    variable_set=coef_s1@Dimnames[[1]][which(coef_s1 != 0)]
    variable_set=variable_set[-1]
    index_elasticnet=as.numeric(sapply(strsplit(variable_set, split='V', fixed=TRUE), function(x) (x[2])))

    num.var=length(index_elasticnet)

    if (num.var==0){
      stop("No variables are selected by ElasticNet, try specify a smaller alpha to relax the penalty")
    }

    G.train.selec=G.train[,index_elasticnet]
    G.test.selec=G.test[,index_elasticnet]

  }

  # generate weights
  weights=rep(0,num.var)

  for(j in 1:J2){
    x=G.train.selec[,j]
    nx=length(x)
    weights[j]=as.numeric((t(y.train)%*%x)/sqrt(nx-1))
  }

  y.test=as.vector(y.test)
  G.test.selec=as.matrix(G.test.selec)

  results=AdapSide(Y=y.test, G=G.test.selec, cor.est=cor.est, lam=lam, weights=weights, pow.param=pow.param)
  return(results)
}

