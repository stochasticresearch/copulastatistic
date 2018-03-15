# Use Mili's method to experiment w/ bivvariate and trivariate

rm(list=ls(all=TRUE))
cat("\014")  

setwd("~/stochasticresearch/copulastatistic/multidimensional_testing")
source('../copulastatistic.R')
source('empcopula.R')
source('conditional_cosdv3d.R')

M = 2000
X1 = rnorm(M)
X2 = rnorm(M)
#X3 = X1 + X2 + rnorm(M,mean=0,sd=0.01)
X3 = X1 + X2

dep1 = cosdv(X1,X2)
dep2 = cosdv(X1,X3)
dep3 = cosdv(X2,X3)

dep4 = cosdv3d(X1,X2,X3)

# compute the copula
zz = empcopula(X1,X2,X3)
cop = unlist(zz[1])
u   = unlist(zz[2])
v   = unlist(zz[3])
w   = unlist(zz[4])

plot(u,cop)
plot(v,cop)
plot(w,cop)

# test computing the conditional CoS
dep1_ccos = conditional_cosdv3d(X1,X2,X3,1,num_bins=10)
dep2_ccos = conditional_cosdv3d(X1,X2,X3,2,num_bins=10)
dep3_ccos = conditional_cosdv3d(X1,X2,X3,3,num_bins=10)

print('With 10 bins!')
print(dep1_ccos)
print(dep2_ccos)
print(dep3_ccos)

# test computing the conditional CoS
dep1_ccos = conditional_cosdv3d(X1,X2,X3,1,num_bins=20)
dep2_ccos = conditional_cosdv3d(X1,X2,X3,2,num_bins=20)
dep3_ccos = conditional_cosdv3d(X1,X2,X3,3,num_bins=20)

print('With 20 bins!')
print(dep1_ccos)
print(dep2_ccos)
print(dep3_ccos)

# test computing the conditional CoS
dep1_ccos = conditional_cosdv3d(X1,X2,X3,1,num_bins=50)
dep2_ccos = conditional_cosdv3d(X1,X2,X3,2,num_bins=50)
dep3_ccos = conditional_cosdv3d(X1,X2,X3,3,num_bins=50)

print('With 50 bins!')
print(dep1_ccos)
print(dep2_ccos)
print(dep3_ccos)

# test computing the conditional CoS
dep1_ccos = conditional_cosdv3d(X1,X2,X3,1,num_bins=100)
dep2_ccos = conditional_cosdv3d(X1,X2,X3,2,num_bins=100)
dep3_ccos = conditional_cosdv3d(X1,X2,X3,3,num_bins=100)

print('With 100 bins!')
print(dep1_ccos)
print(dep2_ccos)
print(dep3_ccos)