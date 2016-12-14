# Code for Figure 7a/b

cat("\014")

n=c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000)
mu=c(0.43,0.27,0.20,0.16,0.14,0.12,0.11,0.09,0.09,0.08,0.07,0.06,0.06,0.06,0.05,0.05,0.05,0.05,0.04,0.04)
sd=c(0.12,0.07,0.05,0.04,0.03,0.03,0.03,0.02,0.02,0.02,0.02,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01)
rhs <- function(x, a,b1)  {
  a*x^b1
} 

xa=n
ya=mu
ds <- data.frame(x = xa, y = ya) 
m <- nls(y ~ rhs(x,a,power), data = ds, start = list(a=1,power = 0), trace = T)
za=sd	
ds <- data.frame(x = xa, y = za) 
m <- nls(y ~ rhs(x,a,power), data = ds, start = list(a=1,power = 0), trace = T)
#for mu curve : approximation : 8.05 n^(-0.74)
#for sd curve : approximation : 2.99 n^(-0.81)
appya= 8.05*n^-0.74
appza=2.99*n^-0.81
par(pty="s")
plot(xa,ya, type = "p", col = "red",xlab="n",ylab=expression(mu),lwd=2,cex=2,pch=16)
lines(xa,appya,lwd=2)
par(pty="s")
plot(xa,za, type = "p", col = "red",xlab="n",ylab=expression(sigma),lwd=2,cex=2,pch=16)
lines(xa,appza,lwd=2)
