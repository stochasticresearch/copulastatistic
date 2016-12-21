source("gentype.r")

n_noise   <- 30
n         <- 500
k1        <- 1
k2        <- 1

power.ace    = power_ace    = array(NA, c(n_types,n_noise))
power.cor    = power_cor    = array(NA, c(n_types,n_noise))
power.dcor   = power_dcor   = array(NA, c(n_types,n_noise))
power.mine   = power_mine   = array(NA, c(n_types,n_noise))
power.rdc    = power_rdc    = array(NA, c(n_types,n_noise))
power.tdc    = power_tdc    = array(NA, c(n_types,n_noise))
power.cos    = power_cos    = array(NA, c(n_types,n_noise))
power.ccor   = power_ccor   = array(NA, c(n_types,n_noise))

for(ttt in c(1,2,3,4,6,7)){
  for(nnn in 1:n_noise){
    fname <- paste("~/ownCloud/PhD/My Papers/CoS/software/results/result", k1, k2, n, ttt, nnn, sep="_")
    #fname <- paste("C:\\Users\\Kiran\\ownCloud\\PhD\\My Papers\\CoS\\software\\results\\result", k1, k2, n, ttt, nnn, sep="_")
    print(paste("Collecting", fname, "..."))
    if(file.exists(fname) == TRUE) load(fname)
    power_ace    [ttt,nnn] = power.ace    [ttt,nnn]
    power_cor    [ttt,nnn] = power.cor    [ttt,nnn]
    power_dcor   [ttt,nnn] = power.dcor   [ttt,nnn]
    power_mine   [ttt,nnn] = power.mine   [ttt,nnn]
    power_rdc    [ttt,nnn] = power.rdc    [ttt,nnn]
    power_tdc    [ttt,nnn] = power.tdc    [ttt,nnn]
    power_cos    [ttt,nnn] = power.cos    [ttt,nnn]
    power_ccor   [ttt,nnn] = power.ccor   [ttt,nnn]
  }
}

power.ace    = power_ace   
power.cor    = power_cor   
power.dcor   = power_dcor  
power.mine   = power_mine
power.rdc    = power_rdc
power.tdc    = power_tdc
power.cos    = power_cos
power.ccor   = power_ccor

library("tikzDevice")

tikz("output_power.tex", standAlone=TRUE, width=7, height=5)

#par(mfrow = c(3,2),lwd=1,mar=rep(0,4),oma=c(4,4,0,0),ps=12)
par(mfrow = c(3,2), lwd=1, oma = c(5,4,0,0) + 0.8,
    mar = c(0,0,1,1) + 0.8, ps=12)

xvals <- c(1:30)/30*3

for(typ in c(1,2,3,4,6,7)) {
  #plot(xvals, power.cor[typ,], pch = 1, col = 1, type = 'b', ylim = c(0,1.1), cex=0.45)
  #points(xvals,power.dcor  [typ,], pch = 2, col = 2, type = 'b')
  plot(xvals, power.dcor[typ,], pch = 2, col = 2, type = 'b', ylim = c(0,1.1))
  points(xvals,power.mine  [typ,], pch = 3, col = 3, type = 'b')
  points(xvals,power.ccor  [typ,], pch = 8, col = 8, type = 'b')
  points(xvals,power.rdc   [typ,], pch = 5, col = 5, type = 'b')
  points(xvals,power.cos   [typ,], pch = 4, col = 4, lwd=3, type = 'b')
  
  if(typ==4)
    legend("topright",c("dCor","TICe", 
                        "Ccor", "RDC", "CoS"), pch =
             c(2,3,8,5,4), col = c(2,3,8,5,4))
  
  #rect(1, .02, 38, .5, col="white")
  rect(0, 0, 0.5, 0.5, col="white")
  xy <- gentype(typ,100,0,0,1)
  xy$x <- (xy$x-min(xy$x))/(max(xy$x)-min(xy$x))*0.5;#+2
  xy$y <- (xy$y-min(xy$y))/(max(xy$y)-min(xy$y))*0.5;#+.05
  ooo  <- order(xy$x)
  points(xy$x[ooo], xy$y[ooo],pch=19,cex=0.3)
}

title(xlab="Noise Level", ylab="Power", outer=TRUE)

dev.off()

system("pdflatex output_power.tex")
