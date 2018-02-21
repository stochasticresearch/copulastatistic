## Script to do some testing of the 3-D CoS

rm(list=ls(all=TRUE))
cat("\014")  

library(MASS)

setwd("~/stochasticresearch/copulastatistic/multidimensional_testing")
source('../copulastatistic.R')

# setup functional dependence tests mentioend by Dr. Mili

M = 500
numMCSim = 100

# Case 1
case1_output = c()
for (mcSimNum in 1:numMCSim) {
  x = runif(M)
  z = x
  y = runif(M)
  case1_output[mcSimNum] = cosdv3d(x,y,z)
}
hist(case1_output)

# Case 2
case2_output = c()
for (mcSimNum in 1:numMCSim) {
  x = runif(M)
  z = (x-0.5)^2
  y = runif(M)
  case2_output[mcSimNum] = cosdv3d(x,y,z)
}
hist(case2_output)

# Case 3
case3_output = c()
sigma = matrix(c(1,0.7,0.7,1),2,2)
mu = c(0,0)
for (mcSimNum in 1:numMCSim) {
  zz= mvrnorm(20,mu,sigma)
  x = zz[,1]
  y = zz[,2]
  z = x
  case3_output[mcSimNum] = cosdv3d(x,y,z)
}
hist(case3_output)
