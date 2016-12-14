# Generate the power curves in Fig. 10

run_reshef <- function(ttt,nnn) {
  
  print(paste("Experiment of type", ttt, "with noise level", nnn))
  
  nbCopRef  = 1
  nbBins    = 10
  nsim      = 500
  n         = 500
  num.noise = 30
  noise     = 3
  val.cor   = val.dcor   = val.mine   = val.ace    = val.rdc   = val.tdc  = val.cos  = val.ccor = rep(NA,nsim)
  val.cor2  = val.dcor2  = val.mine2  = val.ace2   = val.rdc2  = val.tdc2 = val.cos2 = val.ccor2 = rep(NA,nsim)
  power.cor = power.dcor = power.mine = power.ace  = power.rdc = power.tdc = power.cos = power.ccor = array(NA, c(n_types,num.noise))
  
  set.seed(1)
  
  k1 = 1
  k2 = 1
  
  #load cindep
  cindeps = create_cindeps(n,nbBins,nbCopRef)
  #load ctargets
  for(typ in ttt) {
    #if(typ != 5) {
    #  ctargets = create_ctargets(typ,100*n,nbBins,nbCopRef)
    #} else {
    #  ctargets = create_ctargets(1,  100*n,nbBins,nbCopRef)
    #}
    ctargets = create_ctargets(typ,100*n,nbBins,nbCopRef)
    for(l in nnn) {
      for(ii in 1:nsim) {
        xy = gentype(typ, n, l, noise, num.noise)
        x  = runif(n)
        y  = xy$y
        
        val.tdc [ii]   = tdc(x,y,cindeps,ctargets,nbBins)
        val.rdc [ii]   = rdc(x,y)
        val.dcor[ii]   = dcor(x,y)
        val.ace [ii]   = ace(x,y)$rsq
        val.cor [ii]   = abs(cor(x,y))
        val.mine[ii]   = mine(x,y)$TIC
        val.cos[ii]    = cosdv(x,y)
        val.ccor[ii]   = rcd(x,y)
      }
      
      cut.tdc    = quantile(val.tdc, .95)
      cut.rdc    = quantile(val.rdc ,.95)
      cut.dcor   = quantile(val.dcor,.95)
      cut.ace    = quantile(val.ace, .95)
      cut.cor    = quantile(val.cor, .95)
      cut.mine   = quantile(val.mine,.95)
      cut.cos    = quantile(val.cos, .95)
      cut.ccor   = quantile(val.ccor,.95)
      
      for(ii in 1:nsim){
        xy = gentype(typ, n, l, noise, num.noise)
        x  = xy$x
        y  = xy$y
        
        val.tdc2 [ii]   = tdc(x,y,cindeps,ctargets,nbBins)
        val.rdc2 [ii]   = rdc(x,y)
        val.dcor2[ii]   = dcor(x,y)
        val.ace2 [ii]   = ace(x,y)$rsq
        val.cor2 [ii]   = abs(cor(x,y))
        val.mine2[ii]   = mine(x,y)$TIC
        val.cos2[ii]    = cosdv(x,y)
        val.ccor2[ii]   = rcd(x,y)
      }
      
      power.tdc [typ,l]   = sum(val.tdc2  > cut.tdc) /nsim
      power.rdc [typ,l]   = sum(val.rdc2  > cut.rdc) /nsim
      power.dcor[typ,l]   = sum(val.dcor2 > cut.dcor)/nsim
      power.ace [typ,l]   = sum(val.ace2  > cut.ace) /nsim
      power.cor [typ,l]   = sum(val.cor2  > cut.cor) /nsim
      power.mine[typ,l]   = sum(val.mine2 > cut.mine)/nsim
      power.cos [typ,l]   = sum(val.cos2 > cut.cos)/nsim
      power.ccor[typ,l]   = sum(val.ccor2 > cut.ccor)/nsim
      
      print(paste("Type:", typ, "Noise:", l,
                  "ACE:", power.ace [typ,l],
                  "COR:", power.cor [typ,l],
                  "DCO:", power.dcor[typ,l],
                  "TIC:", power.mine[typ,l],
                  "RDC:", power.rdc [typ,l],
                  "TDC:", power.tdc [typ,l],
                  "COS:", power.cos [typ,l],
                  "C2R:", power.ccor[typ,l]))
      cat("\n")
    }
  }
  
  save(power.ace, power.cor, 
       power.dcor, power.mine,
       power.rdc,power.tdc,power.cos,power.ccor,
       file = paste("C:\\Users\\Kiran\\ownCloud\\PhD\\My Papers\\CoS\\software\\results\\result", k1, k2, n, ttt, nnn, sep="_"))
}

#rm(list = ls())
#setwd("~/ownCloud/PhD/My Papers/CoS/software")
setwd("C:\\Users\\Kiran\\ownCloud\\PhD\\My Papers\\CoS\\software")

source("algorithms.R")
source("gentype.r")

#library("doSNOW")
library("foreach")

#cl <- makeCluster(5,type="SOCK")
#registerDoSNOW(cl)

for (ttt in c(1,2,3,4,6,7)) {
  print(paste("typ: ",ttt))
  #foreach(nnn = 1:30) %dopar% {run_reshef(ttt,nnn)}
  foreach(nnn = 1:30) %do% {
    run_reshef(ttt,nnn)
  }
}

#stopCluster(cl)