titles  <- c("Linear", "Quadratic", "Cubic", "Sine", "Cosine",
             "Power 0.25", "Circle", "Step", "Bell", "Two Waves",
             "Logarithm")

n_types <- 11

gentype <- function(typ, n, l, noise, num.noise) {
  x = runif(n)
  
  if(typ==1)  y = x                                                + noise *(l/num.noise)* rnorm(n)
  if(typ==2)  y = 4*(x-.5)^2                                       + noise * (l/num.noise) * rnorm(n)
  if(typ==3)  y = 128*(x-1/3)^3-48*(x-1/3)^3-12*(x-1/3)            + 10* noise  * (l/num.noise) *rnorm(n)
  #if(typ==3)  y = ((-1)^rbinom(n, 1, 0.5))*x                       + noise  * (l/num.noise) *rnorm(n)
  if(typ==4)  y = sin(3*pi*x)                                      + 3*noise * (l/num.noise) *rnorm(n)
  if(typ==5)  y = runif(n) #y = cos(3*pi*x) + 3*noise * (l/num.noise) *rnorm(n)
  if(typ==6)  y = x^(1/4)                                          + noise * (l/num.noise) *rnorm(n)
  if(typ==7)  y = (2*rbinom(n,1,0.5)-1) * (sqrt(1 - (2*x - 1)^2))  + noise/4*l/num.noise *rnorm(n)
  if(typ==8)  y = (x > 0.5)                                        + 7*noise*l/num.noise *rnorm(n)
  if(typ==9)  y = dnorm(x, .5, .1)                                 + 4*noise*l/num.noise*rnorm(n) 
  if(typ==10) y = (2*rbinom(n,1,0.5)-1) * sin(3*pi*x)              + (l/num.noise)*rnorm(n)
  if(typ==11) y = log(x)                                           + 4*noise*l/num.noise*rnorm(n) 
  
  list(x=x,y=y)
}

create_cindeps <- function(n,nbBins=10,nbIndeps) {
  
  cindeps = list()
  for(i in 1:1) {
    # #REF INDEP COPULA
    x <- runif(n)
    y <- runif(n)
    res <- list(x=x,y=y)
    x <- cbind(apply(as.matrix(res$x),2,function(u)rank(u,ties.method="first")/length(u)),1)
    y <- cbind(apply(as.matrix(res$y),2,function(u)rank(u,ties.method="first")/length(u)),1)
    cres = list(x=x,y=y)
    df = data.frame(cres$x[,1],cres$y[,1])
    cindep = hist2d(df,nbins=nbBins,same.scale=TRUE,show=FALSE)
    for(row in 1:nbBins) {
      for(col in 1:nbBins) {
        cindep$counts[row,col] = n/(nbBins*nbBins)
      }
    }
    cindeps[[i]] <- cindep
  }
  cindeps
}

create_ctargets <- function(target,n,nbBins=10,nbTargets) {
  
  ctargets = list()
  for(i in 1:nbTargets) {
    res = gentype(target,n,0,0,30)
    x <- cbind(apply(as.matrix(res$x),2,function(u)rank(u,ties.method="first")/length(u)),1)
    y <- cbind(apply(as.matrix(res$y),2,function(u)rank(u,ties.method="first")/length(u)),1)
    cres = list(x=x,y=y)
    df = data.frame(cres$x[,1],cres$y[,1])
    ctarget = hist2d(df,nbins=nbBins,same.scale=TRUE,show=FALSE)
    ctargets[[i]] <- ctarget
  }
  ctargets
}

