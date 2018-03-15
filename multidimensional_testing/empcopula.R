empcopula <- function(x,y,z) {
  
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
  
  cop
  return(list(cop,u,v,w))
  
}