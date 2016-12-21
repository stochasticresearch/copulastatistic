# Generate Figure 10 in the paper

cat("\014")
source("algorithms.R")

rhoVal = 0.2

i=1
vectx=seq(50,2000,50)
nbrsim=1000
matcos=matrix(nrow=40,ncol=1)
matdcor=matrix(nrow=40,ncol=1)
matmic=matrix(nrow=40,ncol=1)
matrdc=matrix(nrow=40,ncol=1)
matccor=matrix(nrow=40,ncol=1)
nbren=50
while (nbren <=2000) {
  print(nbren)
  s1=0
  s2=0
  s3=0
  s4=0
  s5=0
  for(cptsim in 1:nbrsim) {
    myMvd <- mvdc(copula = normalCopula(rhoVal),
    margins = c("norm", "norm"), paramMargins = list(list(mean = 0,
    sd = 1), list(mean = 0, sd = 1)))
    x.samp <- rMvdc(nbren,myMvd)
    x=x.samp[,1]
    y=x.samp[,2]
    
    s1=s1+cosdv(x,y)
    s2=s2+dcor(x,y)
    s3=s3+mine(x,y)$MICe
    s4=s4+rdc(x,y)
    s5=s5+rcd(x,y)
    #print(cptsim)
  }
  #print(i)
  
  matcos[i,1]=s1/nbrsim
  matdcor[i,1]=s2/nbrsim
  matmic[i,1]=s3/nbrsim
  matrdc[i,1]=s4/nbrsim
  matccor[i,1]=s5/nbrsim
  nbren=nbren+50
  i=i+1
}

cosi=as.vector(matcos)
dcori=as.vector(matdcor)
mati=as.vector(matmic)
rdci=as.vector(matrdc)
ccori=as.vector(matccor)

save.image()

par(pty="s",font=2)
ligne=rep(rhoVal,40)

plot(vectx, dcori, pch = 2, col = 2, type = 'b',xlab="n",ylab="Dependence strength",ylim=c(0.05,1),cex.lab=1.5,cex.axis=1.25)
points(vectx, mati, pch = 3, col = 3, type = 'b')
points(vectx, ccori, pch = 8, col = 8, type = 'b')
points(vectx, rdci, pch = 5, col = 5, type = 'b')
points(vectx, cosi, pch = 4, col = 4, lwd=3, type = 'b')
lines(vectx,ligne,col='black',lwd=2)

legend("topleft",c("dCor","MICe", 
                    "Ccor", "RDC", "CoS"), pch =
         c(2,3,8,5,4), col = c(2,3,8,5,4))
