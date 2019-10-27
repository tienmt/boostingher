source('/data2/thetm/heritability-LM/EigenPrism.R')
library(glmnet)
library(HiLMM)
require(doMC)
registerDoMC(cores=10)

noise = 1

load('/data2/thetm/MA/new.masnps.rda')
masnps.filter = masnps.filter/2
load('/data2/thetm/MA/geneforcausallistMA.RData')
p = ncol(masnps.filter)
n = nrow(masnps.filter)
#names(head(sort(sam.cor,decreasing=T),n))
h2 = 0.2

load('/data2/thetm/MA/ma.cov2.rda')

#c(pbpX,pbp1A,penA)
################################################################
################################################################
Eprism = epr.boost = Enet = MLE = mle.boost = Moment = mm.boost = HERRA = herra.boost = h2aprx= c()
for(ss in 1:30){
spa = 20
spa = sample( 1:p , spa )
beta0 = rep(0,p)
names(beta0) = colnames(masnps.filter)
beta0[spa] = sample(c(-1,-0.5,0.5,1),length(spa),replace = T)
beta0 <- beta0 * sqrt((noise^2)*h2/(1-h2) / as.numeric(crossprod( beta0, ma.cov2) %*% beta0))
y = masnps.filter%*%beta0 + rnorm(n)*noise
h2aprx[ss] = 1 - noise^2/var(y)
#SCREENING
sam.cor <- apply(masnps.filter, 2,function(snp) cor.scale(snp,y) )
xsele = masnps.filter[,sam.cor> quantile(sam.cor, 0.9) ]

#using EigenPrism
Eprism[ss] = EigenPrism(y, xsele, target = 'heritability')$estimate
outlist <- foreach(ii = seq(30)) %dopar%  {
  foreach(i = seq(2), .packages= c("glmnet")) %dopar%  {
    which <- sample(rep(seq(2), length = n)) == i
    xtam = Matrix(xsele[which,], sparse = TRUE)
    tam =  predict(cv.glmnet(xtam, y[which], parallel = T,alpha = .001 ),
                   type = 'nonzero',s = 'lambda.min')[[1]]
    xols1 = xsele[,tam] 
    a = NA
    tryCatch({a <- EigenPrism( y[!which], 
        xols1[!which,], target = 'heritability')$estimate }, error=function(e){ } )
}}
epr.boost[ss]= mean( unlist(outlist) , na.rm=T)

# MLE 
MLE[ss]= estim_herit(y, xsele)$heritability
outlist <- foreach(ii = seq(30)) %dopar%  {
   foreach(i = seq(2), .packages = c("glmnet")) %dopar%  {
     which = sample(rep(seq(2), length = n)) == i
     xols1 = xsele[, predict(cv.glmnet(xsele[which,], y[which],parallel = T,alpha = .001),
                         type = 'nonzero', s = 'lambda.min')[[1]] ]
     a = NA
     tryCatch({a <- her.trunc( estim_herit(y[!which],xols1[!which,])$heritability )  }, 
              error=function(e){ } )
   }}
mle.boost[ss] = mean( unlist(outlist) , na.rm=T)
 
#Dicker moment estimator
Moment[ss] = mm2dicker(y, xsele )
outlist <- foreach(ii = seq(30)) %dopar%  {
   foreach(i = seq(2), .packages = c("glmnet")) %dopar%  {
     which <- sample(rep(seq(2), length = n)) == i
     xols1 = xsele[, predict(cv.glmnet(xsele[which,], y[which], parallel = T,alpha = .01), 
                         type = 'nonzero', s = 'lambda.min')[[1]] ]
     a = NA
     tryCatch({a <- her.trunc(mm2dicker(y[!which], xols1[!which,]))  }, 
              error=function(e){ } )
}}
 mm.boost[ss]= mean( unlist(outlist) , na.rm=T)
 
 # HERRA
 xdat2 = masnps.filter 
 for (itrrs in 1:(log2(p/n)+1) ) {
   model <- cv.glmnet(xdat2 , y , alpha=0, parallel = T)
   coef0 <- abs(predict(model,type="coefficients",s='lambda.min')[-1] )
   xdat2 <- xdat2[, coef0 > quantile(coef0,.5)]
   if(ncol(xdat2) < n ) break
 }
 which = sample(rep(seq(2), length = n)) == 1
 xols1 = xdat2[, predict(cv.glmnet(xdat2[which,], y[which], parallel= T, alpha= 1),
                         type = 'nonzero', s = 'lambda.min')[[1]]  ]
 xols2 = xdat2[, predict(cv.glmnet(xdat2[!which,],y[!which],parallel = T, alpha= 1),
                         type = 'nonzero', s = 'lambda.min')[[1]]  ] 
 a=NA; a <- (summary(lm(y[!which] ~ xols1[!which,] ))$sigma)^2 
 b=NA; b <- (summary(lm(y[which] ~ xols2[which,] ))$sigma)^2
 HERRA[ss]= 1 - mean(na.omit(c(a,b)))/var(y)
 outlist <- foreach(jj = seq(30) ) %dopar%  {
   foreach(i = seq(2), .packages=c("glmnet")) %dopar%  {
     which =  sample(rep(seq(2), length = n)) == i
     xols1 = xdat2[, predict(cv.glmnet(xdat2[which,], y[which], parallel = T,alpha = 1),
                             type = 'nonzero', s = 'lambda.min')[[1]] ]
     a = NA
     tryCatch({a <-(summary(lm(y[!which] ~ xols1[!which,] ))$sigma)^2 
     },  error=function(e){ } )
}}
 herra.boost[ss]= 1 -  mean( unlist(outlist) , na.rm=T)/var(y)
 
 # scaled-lasso
 sig2.lasso = scaled.lasso(xsele , y)^2 
 # ENET
 glmnetcv = cv.glmnet(Matrix(xsele, sparse = T), y, alpha = 0.001, parallel = T)
 tam = predict(glmnetcv, type = 'nonzero', s = 'lambda.min')[[1]] 
 enet = coef( glmnetcv,  s='lambda.min')[-1]
 bSigab = enet[tam] %*% cov2(xsele[,tam]) %*% enet[tam]
 Enet[ss]= bSigab/ (bSigab + sig2.lasso )
 print(ss)
}
save(Eprism,epr.boost, MLE,mle.boost, Enet, Moment,mm.boost, HERRA,herra.boost,h2aprx,
     file = '/data2/thetm/heritability-LM/outsimulations/ma.spa20.h02.noise1.screening90.rdata')


 




 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 which =  sample(rep(seq(2), length = n)) == 1
 tam1 = predict(cv.glmnet(xsele[!which,], y[!which],parallel = T,alpha = .01), 
                type = 'nonzero',s = 'lambda.min')[[1]] 
 tam2 = predict(cv.glmnet(xsele[which,], y[which],parallel = T,alpha = .01),
                type = 'nonzero', s = 'lambda.min')[[1]] 
 x.tam1 = xsele[which,tam1]
 x.tam2 = xsele[!which,tam2]
 
 xtam1.project= x.tam1 %*% solve(crossprod(x.tam1),t(x.tam1))
 y[which] %*% (diag(nrow(xtam1.project))-xtam1.project) %*% y[which]/(n/2 - length(tam1))
 xtam2.project = x.tam2 %*% solve(crossprod(x.tam2),t(x.tam2))
 y[!which] %*% (diag(nrow(xtam1.project))-xtam1.project) %*% y[!which]/(n/2 - length(tam2))
 
 a = mean(xtam1.project.y^2)/(sig2.lasso + mean(xtam1.project.y^2)) 
 xtam2.project.y = x.tam2 %*% solve(crossprod(x.tam2), crossprod(x.tam2, y[!which]))
 b = mean(xtam2.project.y^2)/(sig2.lasso + mean(xtam2.project.y^2)) 
 mean(c(a,b))
 


 


 
 #2-stage strategy
 sam.split = sample(1:nrow(xsele), nrow(xsele)/2)
 tam = predict(cv.glmnet(xsele[sam.split,], 
                         y[sam.split],
                         alpha = 0.01,
                         parallel = T), type = 'nonzero', s = 'lambda.min' )[[1]]
 xnew = xsele[-sam.split,tam]
 S = crossprod(xnew)/n
 her.trunc( mm2dicker(y[-sam.split], xnew )) 
 
 
 