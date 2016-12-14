# Regenerate Figure 8


cat("\014")
source("algorithms.R")
source("copulastatistic.R")
set.seed(1)

vres=matrix(nrow=100,ncol=99)
for ( cpts2 in 1:100) {

  nbri=99
  mycos=rep(0,nbri)
  n=600
  for ( cpts in 1:nbri) {

    print(c(cpts2,cpts))
    myMvd <- mvdc(copula = normalCopula(0.8),
    margins = c("norm", "norm"), paramMargins = list(list(mean = 0,
    sd = 1), list(mean = 0, sd = 1)))
    x.samp <- rMvdc(n,myMvd)

    x=x.samp[,1]
    y=x.samp[,2]
    
    mycos[cpts]=cosdv(x,y)
  }

  mycos=(mycos-mean(mycos))/sd(mycos)
  
  vres[cpts2,]=mycos
  print(cpts2)
}

aaa=rnorm(100)
r=quantile(x,  probs = seq(0.01,0.99,0.01))

vres2=matrix(nrow=100,ncol=99)
for ( aa in 1:100) {
  vres2[aa,]=sort(vres[aa,])
}
vres3=matrix(nrow=99,ncol=100)
for ( aa in 1:99) {
  vres3[aa,]=sort(vres2[,aa])
}

 v1=vres3[,25]
 v2=vres3[,50]
 v3=vres3[,75]
 
 op <- par(mar = c(6,6,4,2) + 0.1, cex.axis=1.5, cex.lab=1.5)
 plot(r,v2,col='blue',typ='p',xlab='theoretical quantiles',ylab='sample quantiles')
 lines(r,v1,col='red',typ='l')
 lines(r,v3,col='green',typ='l')
  grid()
 axis(side = 2, font = 4)
 axis(side = 1, font = 4)
 
#par(pty="s")
#plot(r,v2,col='blue',typ='p',xlab='theoretical quantiles',ylab='sample quantiles')
#lines(r,v1,col='red',typ='l')
#lines(r,v3,col='green',typ='l')
##lines(r,r)
#grid()
