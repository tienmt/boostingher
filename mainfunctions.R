library(compiler)
her.trunc = function(x){
  if(x<0)x=0
  if(x>1)x=1
  return(x)
}
her.trunc <- compiler::cmpfun(her.trunc)
cov2 <- function ( x ) {
  crossprod ( scale (x , TRUE , FALSE ) )/( NROW ( x ) -1)
}
cov2 <- cmpfun(cov2)
cor.scale = function(snp,y){
  a = snp -mean(snp)
  b = y - mean(y)
  abs(a%*%b / sqrt(sum(a^2)*sum(b^2)) )
}
cor.scale <- compiler::cmpfun(cor.scale)

library(cccp)
EigenPrism <- function(y,X,
                       invsqrtSig=NULL,
                       alpha=0.05,
                       target='beta2',
                       zero.ind=c()){
   # Author: Lucas Janson (statweb.stanford.edu/~ljanson)
   # Runs EigenPrism procedure for estimating and generating confidence
   #  intervals for variance components in high-dimensional linear model:
   #       y = X%*%beta + e,   rows of X iid~ N(0,Sig),   e iid~ N(0,sigma^2)
   #  Requires cccp package for solving second order cone optimization.
   #  Note confidence interval endpoints may lie outside parameter domain, so it may be appropriate
   #   to clip them after the fact.
   # 
   # Inputs:
   #  y: response vector of length n (will automatically be centered)
   #  X: n by p design matrix; columns will automatically be centered and scaled to variance 1;
   #      should not contain intercept column, since both y and X will be centered
   #  invsqrtSig: if columns of X not independent, p by p positive definite matrix which is the square-root
   #               of the inverse of Sig, where Sig is the *correlation* matrix of the X (default is identity)
   #  alpha: significance level for confidence interval (default = 0.05)
   #  target: target of estimation/inference
   #		  'beta2' (default) is the squared 2-norm of the coefficient vector: sum(beta^2)
   #           'sigma2' is the noise variance sigma^2
   #           'heritability' is the fraction of variance of y explained by X%*%beta: t(beta)%*%Sig%*%beta/var(y)
   #  zero.ind: vector of which indices of the weight vector w to constrain to zero (default is none)
   #  diagnostics: boolean (default = T) for whether to generate diagnostic plots for the V_i, lambda_i, and w_i
   #  
   # Outputs:
   #  estimate: unbiased estimate of the target (for heritability, only approximately unbiased)
   #  CI: 100*(1-alpha)% confidence interval for target
   
   # Get dimensionality of problem
   n = nrow(X)
   p = ncol(X)
   
   # Transform y and X to proper form
   y = y-mean(y)
   X = scale(X)*n/(n-1)
   if(!is.null(invsqrtSig)) X = X%*%invsqrtSig
   
   # Take singular value decomposition and rescale singular values
   svd = svd(X)
   lambda = svd$d^2/p
   
   # Defined cone-constrained linear problem to optimize weights; [v; w] is vector of optimization variables
   q = c(1,rep(0,n)) #coefficient vector in objective function
   A = rbind(c(0,rep(1,n)),c(0,lambda)) #matrix for linear constraints
   b = c(0,1) #vector for linear constraints
   if(target=='sigma2') b = c(1,0) #switch constraints if target is sigma^2
   # Constrain some weights to be zero if desired
   if(!is.null(zero.ind)){
      A = rbind(A,cbind(rep(0,length(zero.ind)),diag(rep(1,n))[zero.ind,]))
      b = c(b,rep(0,length(zero.ind)))
   }
   # Define second-order cone constraints
   soc1 = socc(diag(c(1/4,rep(1,n))),c(-1/2,rep(0,n)),c(1/4,rep(0,n)),1/2)
   soc2 = socc(diag(c(1/4,lambda)),c(-1/2,rep(0,n)),c(1/4,rep(0,n)),1/2)
   prob = dlp(as.vector(q),A,as.vector(b),list(soc1,soc2))
   
   # Solve optimization problem and extract variables
   opt = cps(prob,ctrl(trace=F))
   v = getx(opt)[1]
   w = getx(opt)[-1]
   
   # Compute estimate and y's variance
   est = sum(w*(t(svd$u)%*%y)^2)
   yvar = sum(y^2)/n
   
   # Compute confidence interval
   CI = est + yvar*sqrt(v)*qnorm(1-alpha/2)*c(-1,1)
   if(target=='heritability'){
      est = est/yvar
      CI = CI/yvar
   }
   # Generate list with results
   result=list()
   result$estimate = est
   result$CI = CI
   #result$para.est = w*(t(svd$u)%*%y)^2
   return(result)
}
scaled.lasso <- function(X, y = NULL, lam0 = NULL) {
  nX = dim(X)[1]
  pX = dim(X)[2]
  if (is.null(lam0)) {
    if (pX > 10^6) {
      lam0 = "univ"
    } else lam0 = "quantile"
  }
  if (lam0 == "univ" | lam0 == "universal")lam0 = sqrt(2 * log(pX)/nX)
  if (lam0 == "quantile") {
    L = 0.1
    Lold = 0
    while (abs(L - Lold) > 0.001) {
      Lold = L
      L = -qnorm(min((L^4 + 2 * L^2)/pX , 0.99))
      L = (L + Lold)/2
    }
    if (pX == 1) 
      L = 0.5
    lam0 = sqrt(2/nX) * L
  }
  
  sigmaint = sqrt(sum(y^2)/nX)
  flag = 0
  sigmanew = sigmaint + 5
  while (abs(sigmaint - sigmanew) > 1e-04 & flag <= 100) {
    flag = flag + 1
    sigmaint = sigmanew
    lam = lam0 * sigmaint
    hbeta = coef(glmnet::glmnet(X,y,intercept = F, 
                                standardize = F,lambda = lam))[-1]
    sigmanew = sqrt(sum((y-X%*%hbeta)^2)/nX)
  }
  as.numeric(sigmanew)
}
mm2dicker = function(y,x,CI = FALSE){
  n = nrow(x)
  p = ncol(x)
  S = crossprod(x)/n
  trcxtx.n = sum(diag(S))
  m1 = trcxtx.n/p
  m2 = sum(S^2)/p  - trcxtx.n^2/n/p
  sum.Xty2 =  sum((t(x)%*%y)^2)
  tau2.dicker = sum.Xty2*m1/( n*(n+1)*m2 ) -  mean(y^2)*p*m1^2 / ((n+1)*m2)
  sig2.dicker = (1+ p*m1^2/m2/(n+1))*mean(y^2) - sum.Xty2*m1/( n*(n+1)*m2 )
  if(CI ==FALSE) return( her.trunc(tau2.dicker/ mean(y^2)) )
  if(CI==TRUE){
  m3 = sum(diag(S^3))/p  +  2*trcxtx.n^3/(p*n^2) - 3*trcxtx.n*sum(diag(S^2))/p/n
  phi0 = 2*( p*m1^2*mean(y^2)^2/n/m2  + 2*m1*m3*(sig2.dicker*tau2.dicker + tau2.dicker^2)/m2^2 -  tau2.dicker^2) / mean(y^2)^2
  confident = tau2.dicker/ mean(y^2) + phi0*qnorm(1-0.05/2)*c(-1,1)/sqrt(n) 
  list(est = her.trunc(tau2.dicker/ mean(y^2) ),CI = confident)
  }
}

EigenPrism <- cmpfun(EigenPrism)
scaled.lasso <- cmpfun(scaled.lasso)
mm2dicker <- cmpfun(mm2dicker)
