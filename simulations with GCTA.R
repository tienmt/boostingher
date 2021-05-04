#corrplot::corrplot(cov2(x[,10000:10200]),method = "shade",tl.pos='n',order = "hclust",col = gray.colors(100))
cov2 <- function ( x ) {
  crossprod ( scale (x , TRUE , FALSE ) )/( NROW ( x ) -1)
}
cov2 <- compiler::cmpfun(cov2)
cor.scale = function(snp,y){
  a = snp -mean(snp)
  b = y - mean(y)
  abs(a%*%b / sqrt(sum(a^2)*sum(b^2)) )
}
cor.scale <- compiler::cmpfun(cor.scale)


library(glmnet)
require(doMC)
registerDoMC(cores=10)

load('/data2/thetm/MA/new.masnps.rda')
masnps.filter = masnps.filter/2
#
p = ncol(masnps.filter)
n = nrow(masnps.filter)


#load('/data2/thetm/MA/geneforcausallistMA.RData')
noise = 1
h2 = 0.5
spa = sample(p,100,replace = F)  #sample(c(pbpX,pbp1A,penA),100)#
beta0 = rep(0,p)
names(beta0) = colnames(masnps.filter)
beta0[spa] = rnorm(length(spa))
beta0 <- beta0 * sqrt((noise^2)*h2/(1-h2) / as.numeric(beta0[spa]%*%cov2(masnps.filter[,spa])%*%beta0[spa]) )



################################################################
################################################################
Enet = HERRA = herra.boost = h2aprx= c()
for(ss in 1:30){
  y = gcta5000phen[,ss]
  #SCREENING
  sam.cor <- apply(masnps.filter, 2,function(snp) cor.scale(snp,y) )
  xsele = masnps.filter[,sam.cor> quantile(sam.cor, 0.60) ]
  p = ncol(xsele)
  # HERRA
  xdat2 = xsele 
  for (itrrs in 1:(log2(p/n)+1) ) {
    model <- cv.glmnet(xdat2 , y , alpha=0, parallel = T)
    coef0 <- abs(predict(model,type="coefficients",s='lambda.min')[-1] )
    xdat2 <- xdat2[, coef0 > quantile(coef0,.5)]
    if(ncol(xdat2) < n ) break
  }
  which <- sample(rep(seq(2), length = n)) == 1
  xols1 = xdat2[, predict(cv.glmnet(xdat2[which,], y[which], alpha= 1),
                          type = 'nonzero', s = 'lambda.min')[[1]]  ]
  xols2 = xdat2[, predict(cv.glmnet(xdat2[!which,],y[!which], alpha= 1),
                          type = 'nonzero', s = 'lambda.min')[[1]]  ] 
  a=NA; a <- (summary(lm(y[!which] ~ xols1[!which,] ))$sigma)^2 
  b=NA; b <- (summary(lm(y[which] ~ xols2[which,] ))$sigma)^2
  HERRA[ss]= 1 - mean(na.omit(c(a,b)))/var(y)
  #boosting HERRA
  outlist <- foreach(jj = seq(50) ) %dopar%  {
    foreach(i = seq(2), .packages=c("glmnet")) %dopar%  {
      which <- sample(rep(seq(2), length = n)) == i
      xols1 = xdat2[, predict(cv.glmnet(xdat2[which,], y[which], parallel = T,alpha = 1),
                              type = 'nonzero', s = 'lambda.min')[[1]] ]
      a = NA
      tryCatch({a <-(summary(lm(y[!which] ~ xols1[!which,] ))$sigma)^2 
      },  error=function(e){ } )
    }}
  herra.boost[ss]= 1 - mean( unlist(outlist), na.rm=T)/var(y)
  
  # ENET
  glmnetcv = cv.glmnet(xsele, y, alpha = 0.01,parallel = T)
  tam = predict(glmnetcv, type = 'nonzero', s = 'lambda.min')[[1]] 
  enet = coef( glmnetcv,  s='lambda.min')[-1]
  Enet[ss]= enet[tam] %*% cov2(xsele[,tam]) %*% enet[tam]/ var(y)
  print(ss)
}

save(Enet, HERRA,herra.boost,h2aprx,
     file = '/data2/thetm/heritability-LM/outsimulations/revised.gcta5000h05.screen60.rda')



mean(Enet);mean(HERRA);mean(herra.boost)
sd(Enet);sd(HERRA);sd(herra.boost)

Enet = c()
for(ss in 1:30){
  y = masnps.filter%*%beta0 + rnorm(n)*noise
  sam.cor <- apply(masnps.filter, 2,function(snp) cor.scale(snp,y) )
  xsele = masnps.filter[,sam.cor> quantile(sam.cor, 0.60) ]
  
  glmnetcv = cv.glmnet(xsele, y, alpha = 0.01, parallel = T)
  tam = predict(glmnetcv, type = 'nonzero', s = 'lambda.min')[[1]] 
  enet = coef( glmnetcv,  s='lambda.min')[-1]
  Enet[ss] = enet[tam] %*% cov2(xsele[,tam]) %*% enet[tam]/ var(y)
  print(ss)
}


