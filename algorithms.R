# Required Packages

# The following packages can be installed from CRAN
library(acepack)
library(minerva)
library(kernlab)
library(energy)
library(gplots)
library(emdist)
library(copula)

# The following package should be installed from the link:
# https://github.com/liyi-1989/rcd
library(rcd)

printf <- function(...) invisible(print(sprintf(...)))

lcg <- function(a,c,m,n) {
  ThePts = matrix(nrow =2 * n,ncol=1)
  ThePts[1,1] = 0
  for (nn in  1 : (2 * n - 1))
  {
    ThePts[nn+1,1] = ((a * ThePts[nn,1])+c)%% m
    
  }
  
  ThePts = ThePts / m;
  
  x = matrix(nrow = n,ncol=1)
  y = matrix(nrow = n,ncol=1)
  for (nn in  1: n)
  {
    indx = 2 * nn - 1
    x[nn, 1] = sqrt(-2 * log(ThePts[indx,1])) * cos(2 * pi * ThePts[indx+1,1])
    y[nn, 1] = sqrt(-2 * log(ThePts[indx,1])) * sin(2 * pi * ThePts[indx+1,1])
  }
  x[1,1] = 0; y[1,1] = 0;   # log(0) is infinity, override.
  df <- data.frame(x, y)
  return(df)
}

# Code from Lopez-Paz
# https://github.com/lopezpaz/randomized_dependence_coefficient
computeKernelMatrix <- function(sample) {
  n        <- nrow(sample)
  Q        <- matrix(apply(sample^2, 1, sum), n, n)
  distance <- Q + t(Q) - 2 * sample %*% t(sample)
  exp(-sigest(sample,scale=NULL)[2]*distance)
}

# Code from Lopez-Paz
# https://github.com/lopezpaz/randomized_dependence_coefficient
hsic <- function(sampleX, sampleY) {
  N  <- nrow(as.matrix(sampleX))
  K  <- computeKernelMatrix(as.matrix(sampleX))
  L  <- computeKernelMatrix(as.matrix(sampleY))
  KH <- K - 1 / N * matrix(apply(K, 2, sum), N, N)
  LH <- L - 1 / N * matrix(apply(L, 2, sum), N, N)
  1 / N * sum(sum(KH * t(LH)))
}

# Code from Lopez-Paz
# https://github.com/lopezpaz/randomized_dependence_coefficient
hsiccop <- function(sampleX, sampleY) {
  sampleX <- apply(as.matrix(sampleX),2,function(u) ecdf(u)(u))
  sampleY <- apply(as.matrix(sampleY),2,function(u) ecdf(u)(u))
  hsic(sampleX,sampleY)
}

# Code from Lopez-Paz
# https://github.com/lopezpaz/randomized_dependence_coefficient
rdc <- function(x,y,k=20,s=1/6,f=sin) {
  x <- cbind(apply(as.matrix(x),2,function(u)rank(u)/length(u)),1)
  y <- cbind(apply(as.matrix(y),2,function(u)rank(u)/length(u)),1)
  x <- s/ncol(x)*x%*%matrix(rnorm(ncol(x)*k),ncol(x))
  y <- s/ncol(y)*y%*%matrix(rnorm(ncol(y)*k),ncol(y))
  cancor(cbind(f(x),1),cbind(f(y),1))$cor[1]
}

# code from Data-grapple
# https://www.datagrapple.com/Tech/optimal-copula-transport.html
tdc <- function(x,y,cindeps,ctargets,nbBins=10) {
  #transform data to empirical copula cdata
  x <- cbind(apply(as.matrix(x),2,function(u)rank(u,ties.method="first")/length(u)),1)
  y <- cbind(apply(as.matrix(y),2,function(u)rank(u,ties.method="first")/length(u)),1)  
  cres = list(x=x,y=y)
  df_cdata = data.frame(cres$x[,1],cres$y[,1])
  cdata = hist2d(df_cdata,nbins=nbBins,same.scale=TRUE,show=FALSE)
  
  
  #compute distance between cindep and cdata
  dists = list()
  for(i in 1:length(cindeps)) {
    cindep <- cindeps[[i]]
    dists[[i]] = emd2d(cindep[1]$counts/length(x),cdata[1]$counts/length(x),xdist=1/nbBins,ydist=1/nbBins,dist="manhattan")
  }
  indep2data = min(unlist(dists))
  
  #compute min distance between cdata and ctargets
  dists = list()
  for(i in 1:length(ctargets)) {
    ctarget <- ctargets[[i]]
    dists[[i]] = emd2d(cdata[1]$counts/length(x),ctarget[1]$counts/length(x),xdist=1/nbBins,ydist=1/nbBins,dist="manhattan")
  }
  data2dep = mean(unlist(dists))
  
  #compute TDC
  (indep2data / (indep2data + data2dep))
}
