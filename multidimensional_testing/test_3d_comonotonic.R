## Script to do some testing of the 3-D CoS

rm(list=ls(all=TRUE))
cat("\014")  

setwd("~/stochasticresearch/copulastatistic/multidimensional_testing")
source('../copulastatistic.R')

M = 500
x1 = runif(M)
x2 = runif(M)
x3 = 1 - (x1 + x2)

zz = cosdv3d(x1,x2,x3)

xx = c(0,0.005,0.01,0.015,0.03)
yy = c(0.07,0.065,0.06,0.055,0.045)
zz = c(0.93,0.93,0.93,0.93,0.925)

dep1 = cosdv(xx,yy)
dep2 = cosdv(xx,zz)
dep3 = cosdv(yy,zz)