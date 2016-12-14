# Regenerate table-2

cat("\014")
source("algorithms.R")
source("copulastatistic.R")

set.seed(123)

ss = c(100,500,1000,2000,3000)
nsim = 1000

alpha = .01
I = c(alpha/2, 1-alpha/2)

# TODO: fill these vectors in w/ results from TABLE1's run!
mu0 = c(0.28,0.08,0.04,0.02,0.02)
std0 = c(0.08,0.02,0.01,0.01,0.01)

idx = 1;
for(M in ss) {
  
  ga1Results = 1:nsim
  ga2Results = 1:nsim
  
  for(ii in 1:nsim) {
    ga1Data = rCopula(M, normalCopula(0.1))
    ga2Data = rCopula(M, normalCopula(0.3))
    
    ga1Results[ii] = cosdv(ga1Data[,1], ga1Data[,2])
    ga2Results[ii] = cosdv(ga2Data[,1], ga2Data[,2])
  }
  
  # get the null's mu0 and std0
  mu0Val = mu0[idx]
  std0Val = std0[idx]
  
  # compute the actual mean
  muGa1 = mean(ga1Results)
  muGa2 = mean(ga2Results)
  
  q = mu0Val + qt(I, df=10000) * std0Val
  # compute the Type II Error
  p1 = pt((q - muGa1)/std0Val, df=10000)
  type2ErrGa1 = diff(p1)
  
  p2 = pt((q - muGa2)/std0Val, df=10000)
  type2ErrGa2 = diff(p2)
  
  printf("*****************************")
  printf("M = %d Gaussian(0.1).Type2Err=%0.02f Gaussian(0.3).Type2Err=%0.02f", M, type2ErrGa1, type2ErrGa2)
  printf("*****************************")
  
  idx = idx + 1;
}