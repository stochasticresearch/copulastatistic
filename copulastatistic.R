# Algorithm for computing the Copula Statistic (CoS)
# See paper "A Copula Statistic for Measuring Nonlinear Multivariate Dependence"
# for more information on the algorithm

# *** NOTE ***
# Code written by Mohsen Ben Hassine (mohsenmbh851@gmail.com).  Feel free to 
# correspond directly with the author, or myself (Kiran Karra - kiran.karra@gmail.com)
# if you have any issues with the CoS implementations provided below!

# *** TODO ***
# [ ] - Generalize cosdv, cosdv3d, cosdv4d such that the cosdv function
#       can detect the dimensionality of the and compute all pseudo-observations
#       without having the user call the appropriate function

cosdv <- function(x, y) {
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
  res
}

cosdv3d<- function(x, y,z) {
  n=length(x)
  vect <- order(x , y,z )
  x1=x
  y1=y
  z1=z
  for ( i in 1:n) {
    x1[i]=x[vect[i]]
    y1[i]=y[vect[i]]
    z1[i]=z[vect[i]]
  }
  x=x1
  y=y1
  z=z1
  
  data<-list(u=x,v=y,w=z)
  pseudodata <- sapply(data, rank, ties.method = "random") / (n+1 )
  u=pseudodata[,1]
  v=pseudodata[,2]
  w=pseudodata[,3]
  
  cop=u
  for(i in 1:n) {
    cpt=0
    for(j in 1:i) {
      if (( u[i]>=u[j]) && ( v[i]>=v[j]) && ( w[i]>=w[j])) {
        cpt=cpt+1
      }
    }
    cop[i]=cpt/(n + 1)
  }
  min=pmin(u,v,w)
  max=pmax((u+v+w)-2,0)
  pip=u*v*w
  
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
          if (cop[i]<pip[i]) {
            w=(pip[i]-cop[i])/(pip[i]-max[i])
          }
          else {
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
  res
}

cosdv4d<- function(x,y,z,t) {
  n=length(x)
  vect <- order(x,y,z,t )
  x1=x
  y1=y
  z1=z
  t1=t
  
  for ( i in 1:n) {
    x1[i]=x[vect[i]]
    y1[i]=y[vect[i]]
    z1[i]=z[vect[i]]
    t1[i]=t[vect[i]]
  }
  x=x1
  y=y1
  z=z1
  t=t1
  data<-list(u=x,v=y,w=z,s4=t)
  pseudodata <- sapply(data, rank, ties.method = "random") / (n+1 )
  u=pseudodata[,1]
  v=pseudodata[,2]
  w=pseudodata[,3]
  s4=pseudodata[,4]
  cop=u
  for(i in 1:n) {
    cpt=0
    for(j in 1:i) {
      if (( u[i]>=u[j]) && ( v[i]>=v[j]) && ( w[i]>=w[j]) && ( s4[i]>=s4[j]) ) {
        cpt=cpt+1
      }
    }
    cop[i]=cpt/(n + 1)
  }
  min=pmin(u,v,w,s4)
  max=pmax((u+v+w+s4)-3,0)
  pip=u*v*w*s4
  
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
      if (condd && indep) w=0
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
          if (cop[i]<pip[i]) {
            w=(pip[i]-cop[i])/(pip[i]-max[i])
          }
          else {
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
  res
}