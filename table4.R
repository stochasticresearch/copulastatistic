# Regenerate the results for Table-4

cat("\014")
source("algorithms.R")
source("copulastatistic.R")

set.seed(123)

ss = c(100,500,1000,2000)
nsim = 1000;

for(M in ss) {
  
  ga1Results = 1:nsim
  ga2Results = 1:nsim
  gu1Results = 1:nsim
  gu2Results = 1:nsim
  cl1Results = 1:nsim
  cl2Results = 1:nsim
  
  for(ii in 1:nsim) {
    ga1Data = rCopula(M, normalCopula(0.1))
    ga2Data = rCopula(M, normalCopula(0.3))
    
    gu1Data <- rCopula(M, gumbelCopula(1.08))
    gu2Data <- rCopula(M, gumbelCopula(1.26))
    
    cl1Data <- rCopula(M, claytonCopula(0.15))
    cl2Data <- rCopula(M, claytonCopula(0.51))
    
    ga1Results[ii] = cosdv(ga1Data[,1], ga1Data[,2])
    ga2Results[ii] = cosdv(ga2Data[,1], ga2Data[,2])
    
    gu1Results[ii] = cosdv(gu1Data[,1], gu1Data[,2])
    gu2Results[ii] = cosdv(gu2Data[,1], gu2Data[,2])
    
    cl1Results[ii] = cosdv(cl1Data[,1], cl1Data[,2])
    cl2Results[ii] = cosdv(cl2Data[,1], cl2Data[,2])
    
  }
  
  printf('*********************\n')
  
  printf("M=%d Gauss(0.1).mu=%0.02f Gauss(0.1).std=%0.02f", M, mean(ga1Results), sd(ga1Results))
  printf("M=%d Gauss(0.3).mu=%0.02f Gauss(0.3).std=%0.02f", M, mean(ga2Results), sd(ga2Results))
  
  printf("M=%d Gumbel(1.08).mu=%0.02f Gumbel(1.08).std=%0.02f", M, mean(gu1Results), sd(gu1Results))
  printf("M=%d Gumbel(1.26).mu=%0.02f Gumbel(1.26).std=%0.02f", M, mean(gu2Results), sd(gu2Results))
  
  printf("M=%d Clayton(0.15).mu=%0.02f Clayton(0.15).std=%0.02f", M, mean(cl1Results), sd(cl1Results))
  printf("M=%d Clayton(0.51).mu=%0.02f Clayton(0.51).std=%0.02f", M, mean(cl2Results), sd(cl2Results))
  
  printf('*********************\n')
}
