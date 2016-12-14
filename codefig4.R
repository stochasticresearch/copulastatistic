# Generate Figure 4 for the Paper

library(copula)
cat("\014")
set.seed(1)
myMvd <- mvdc(copula = normalCopula(0.5),
margins = c("norm", "norm"), paramMargins = list(list(mean = 0,
sd = 1), list(mean = 0, sd = 1)))
x.samp <- rMvdc(100,myMvd)
x=x.samp[,1]
y=x.samp[,2]

n=length(x)
vect <- order(x , y )
x1=x
y1=y
for ( i in 1:n) {
  x1[i]=x[vect[i]]
  y1[i]=y[vect[i]]
}
x=x1
y=y1
data<-list(u=x,v=y)
pseudodata <- sapply(data, rank, ties.method = "random") / (n+1 )
u=pseudodata[,1]
v=pseudodata[,2]
cop=u
for(i in 1:n) {
  cpt=0
  for(j in 1:i) {
    if (( u[i]>=u[j]) && ( v[i]>=v[j])) {
      cpt=cpt+1
    }
  }
  cop[i]=cpt/(n + 1)
}
min=pmin(u,v)
max=pmax((u+v)-1,0)
pip=u*v

score=0
wanc=1

ss=cop[1]
if (cop[1]==pip[1]) {
  wanc=0
}
if ((cop[1]==min[1])||(cop[1]==max[1])||(cop[1]==cop[2])) {
  wanc=1
}

indep=(mean(abs(diff(cop)))>0.12)

for ( i in 2:n) {
  ss=c(ss,cop[i])
  if (( is.unsorted (ss))&&( is.unsorted (rev(ss)))) {
    cond1=(cop[i-1]==cop[i-2]) 
    if (cond1) {
      j=i-2
    }
    else {
      j=i-1  
    }
    if (cop[j]>pip[j]) {
      w=(cop[j]-pip[j])/(min[j]-pip[j])
    }
    else {
      if (cop[j]<pip[j]) {
        w=(pip[j]-cop[j])/(pip[j]-max[j])
      }
      else {
        w=0
      }
    }
    if (((round(cop[i-1],2)==round(min[i-1],2))||(round(cop[i-1],2)==round(max[i-1],2))|| (cop[i-1]==(1/(n+1))))&& (length(ss)>4 )) {
      w=1
    }
    if (((cop[i]==cop[i-2])||(cop[i-1]==cop[i-2]))&& (length(ss)>4 )) {
      w=1
    }
    if (((round(pip[i-1],2)==round(min[i-1],2))||(round(pip[i-1],2)==round(max[i-1],2)))&& (length(ss)>4 )) {
      w=1
    }
    condd=(round(cop[i-1],2)==round(pip[i-1],2))
    if (condd && indep) {
      w=0
    }
    score=score+(length(ss)-1)*(w+wanc)/2
    wanc=w
    ss=cop[i]
  }
  else {
    if (i==n) {
      if (cop[i]>pip[i]) {
        w=(cop[i]-pip[i])/(min[i]-pip[i])
      }
      else {
        if (cop[i]<pip[i]){
          w=(pip[i]-cop[i])/(pip[i]-max[i])
        }
        else{
          w=0
        }
      }
      if (((cop[i-1]==min[i-1])||(cop[i-1]==max[i-1]))&& (length(ss)>=4 )) {
        w=1
      }
      score=score+(length(ss)-1)*(w+wanc)/2
    }
  }
}
res<-(score+1)/n

op <- par(mar = c(6,6,4,2) + 0.1, cex.axis=1.5, cex.lab=2)
plot(u,cop, type = "p", col = "blue",xlab="u",ylab=expression('C'['uv']),xlim=c(0,0.3),ylim=c(0,0.3))
lines(u,pi, type = "l", col = "red")
lines(u,min, type = "l",lty=3, col = "black")
lines(u,max, type = "l", lty=2,col = "green")
axis(side = 2, font = 4)
axis(side = 1, font = 4)
#legend(3, 0.9, c("Cuv", "PI", "M","W"), col = c("blue", "red","black", "green"),
#       text.col = "black", lty = c(1,1,3,2),
#       merge = TRUE, bg = "gray90", cex=1.5)
#legend("topleft", c(expression('C'['uv']),expression(paste(Pi)),"M","W"), col = c("blue", "red","black", "green"),
#       lty = c(1,1,3,2), cex=2.5)

