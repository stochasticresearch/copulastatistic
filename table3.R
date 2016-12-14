# Regenerate the results for Table-4

cat("\014")
source("algorithms.R")
source("copulastatistic.R")

set.seed(123)

ss = c(100,500,1000,2000,3000,5000)
nsim = 1000;

for(M in ss) {
  
  res1 = 1:nsim
  res2 = 1:nsim
  res3 = 1:nsim
  
  for(ii in 1:nsim) {
    x1 = rnorm(M)
    y1 = sin(x1)
    
    x2 = rnorm(M)
    y2 = sin(5*x2)
    
    x3 = rnorm(M)
    y3 = sin(14*x3)
    
    res1[ii] = cosdv(x1,y1)
    res2[ii] = cosdv(x2,y2)
    res3[ii] = cosdv(x3,y3)
  }
  
  printf('*********************')
  printf("Dep1 --> M=%d mu=%0.02f std=%0.02f", M, mean(res1), sd(res1))
  printf("Dep2 --> M=%d mu=%0.02f std=%0.02f", M, mean(res2), sd(res2))
  printf("Dep3 --> M=%d mu=%0.02f std=%0.02f", M, mean(res3), sd(res3))
  printf('*********************')
}
