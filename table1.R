# Regenerate the results for Table-1

cat("\014")
source("algorithms.R")
source("copulastatistic.R")

set.seed(123)

ss = c(100,500,1000,2000,3000)
nsim = 1000;

for(M in ss) {
  
  ga1Results = 1:nsim
  gu1Results = 1:nsim
  cl1Results = 1:nsim
  
  for(ii in 1:nsim) {
    ga1Data = rCopula(M, normalCopula(0))
    gu1Data <- rCopula(M, gumbelCopula(1.000000000000001))  # we do this to suppress warning messages easily :0
    cl1Data <- rCopula(M, claytonCopula(0.0000000000000000000001)) # we do this to suppress warning messages easily :0
    
    ga1Results[ii] = cosdv(ga1Data[,1], ga1Data[,2])
    gu1Results[ii] = cosdv(gu1Data[,1], gu1Data[,2])
    cl1Results[ii] = cosdv(cl1Data[,1], cl1Data[,2])
  }
  
  printf('*********************\n')
  
  printf("M=%d Gauss(0).mu=%0.02f Gauss(0).std=%0.02f", M, mean(ga1Results), sd(ga1Results))
  printf("M=%d Gumbel(1).mu=%0.02f Gumbel(1).std=%0.02f", M, mean(gu1Results), sd(gu1Results))
  printf("M=%d Clayton(0).mu=%0.02f Clayton(0).std=%0.02f", M, mean(cl1Results), sd(cl1Results))

  printf('*********************\n')
}
