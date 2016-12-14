# Generate the Equitability Curves in Figure 9

cat("\014")
source("algorithms.R")
source("copulastatistic.R")

set.seed(1)
xi=seq(0.01,0.99,0.01)

resf1=matrix(nrow=99,ncol=100)
resf2=matrix(nrow=99,ncol=100)
resf3=matrix(nrow=99,ncol=100)
resf4=matrix(nrow=99,ncol=100)
resf5=matrix(nrow=99,ncol=100)
resf6=matrix(nrow=99,ncol=100)
resf7=matrix(nrow=99,ncol=100)
resf8=matrix(nrow=99,ncol=100)
resf9=matrix(nrow=99,ncol=100)
resf10=matrix(nrow=99,ncol=100)
R=0.01
Rinv=1/R
samplesize=250

for (i in 1:99) {
  for (j in 1:100) {
    x=runif(samplesize)
    y=cos(14*pi*x)
    noise=rnorm(samplesize,0,sqrt(var(y)*(Rinv-1)))
    y=y+noise
    resf1[i,j]=cosdv(x,y)
    
    x=runif(samplesize)
    y=sin(5*pi*x*(1 + x))
    noise=rnorm(samplesize,0,sqrt(var(y)*(Rinv-1)))
    y=y+noise
    resf2[i,j]=cosdv(x,y)
    
    x=runif(samplesize,-1.3,1.1)
    y=41*(4*x^3 + x^2 - 4*x)
    noise=rnorm(samplesize,0,sqrt(var(y)*(Rinv-1)))
    y=y+noise
    resf3[i,j]=cosdv(x,y)
    
    x=runif(samplesize,0,10)
    y=2^x
    noise=rnorm(samplesize,0,sqrt(var(y)*(Rinv-1)))
    y=y+noise
    resf4[i,j]=cosdv(x,y)
    
    x=runif(samplesize)
    y=x
    noise=rnorm(samplesize,0,sqrt(var(y)*(Rinv-1)))
    y=y+noise
    resf5[i,j]=cosdv(x,y)
    
    x=runif(samplesize)
    y=10*sin(10.6*(2*x - 1)) + 11/10 *(2*x - 1)
    noise=rnorm(samplesize,0,sqrt(var(y)*(Rinv-1)))
    y=y+noise
    resf6[i,j]=cosdv(x,y)
    
    x=runif(samplesize)
    y=sin(10*pi*x) + x
    noise=rnorm(samplesize,0,sqrt(var(y)*(Rinv-1)))
    y=y+noise
    resf7[i,j]=cosdv(x,y)
    
    x=runif(samplesize,-0.5,0.5)
    y=4*(x^2)
    noise=rnorm(samplesize,0,sqrt(var(y)*(Rinv-1)))
    y=y+noise
    
    resf8[i,j]=cosdv(x,y)
    
    x=runif(samplesize)
    y=sin(16*pi*x)
    noise=rnorm(samplesize,0,sqrt(var(y)*(Rinv-1)))
    y=y+noise
    resf9[i,j]=cosdv(x,y)
    
    x=runif(samplesize)
    y=sin(6*pi*x*(1 + x))
    noise=rnorm(samplesize,0,sqrt(var(y)*(Rinv-1)))
    y=y+noise
    resf10[i,j]=cosdv(x,y)
  }
  
  R=R+0.01
  Rinv=1/R
  
  print(i)
}

j=1
par(pty="s")
plot(xi,resf1[,j], type = "p", col = 1,xlab='R2',ylab='CoS',ylim=c(0,1))
lines(xi,resf2[,j], type = "p", col = 8)
lines(xi,resf3[,j], type = "p", col = 3)
lines(xi,resf4[,j], type = "p", col = 4)
lines(xi,resf5[,j], type = "p", col = 5)
lines(xi,resf6[,j], type = "p", col = 6)
lines(xi,resf7[,j], type = "p", col = 7)
lines(xi,resf8[,j], type = "p", col = 2)
lines(xi,resf9[,j], type = "p", col = "#ccFF00")
lines(xi,resf10[,j], type = "p", col = "#FD3F92")


for(j in 2:99) {
  lines(xi,resf1[,j], type = "p", col = 1)
  lines(xi,resf2[,j], type = "p", col = 8)
  lines(xi,resf3[,j], type = "p", col = 3)
  lines(xi,resf4[,j], type = "p", col = 4)
  lines(xi,resf5[,j], type = "p", col = 5)
  lines(xi,resf6[,j], type = "p", col = 6)
  lines(xi,resf7[,j], type = "p", col = 7)
  lines(xi,resf8[,j], type = "p", col = 2)
  lines(xi,resf9[,j], type = "p", col = "#ccFF00")
  lines(xi,resf10[,j], type = "p", col = "#FD3F92")
}



