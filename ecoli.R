# Generate all the required data for the ecoli expression processing

cat("\014")
source("algorithms.R")
source("copulastatistic.R")

library(minet)
library(readr)

cosPairwiseFilename = "./ecoli_results/cos.csv"
micePairwiseFilename = "./ecoli_results/mice.csv"
dcorPairwiseFilename = "./ecoli_results/dcor.csv"
rdcPairwiseFilename = "./ecoli_results/rdc.csv"
ccorPairwiseFilename = "./ecoli_results/ccor.csv"

# read in the data
ecoliData = read_csv("./ecoli_results/data.csv")
numFeatures = ncol(ecoliData)

# compute the pairwise dependencies for all 5 metrics 
# CoS, MICe, dCor, RDC, cCor
# if the processing hasn't already been done
if(!file.exists(cosPairwiseFilename)) {
  printf("Computing CoS Pairwise Dependencies for eColi Data")
  pairwiseMat = matrix(0, nrow=numFeatures, ncol=numFeatures)
  for(ii in c(1:numFeatures)) {
    for(jj in c(ii:numFeatures)) {
      if(ii!=jj) {
        # a symetric matrix ...
        pairwiseMat[ii,jj] = cosdv(ecoliData[[ii]], ecoliData[[jj]])
        pairwiseMat[jj,ii] = cosdv(ecoliData[[ii]], ecoliData[[jj]])
      }
    }
  }
  # write to file
  write.table(pairwiseMat, file = cosPairwiseFilename, append = FALSE, row.names = F, sep = ",");
}
if(!file.exists(micePairwiseFilename)) {
  printf("Computing MICe Pairwise Dependencies for eColi Data")
  pairwiseMat = matrix(0, nrow=numFeatures, ncol=numFeatures)
  for(ii in c(1:numFeatures)) {
    for(jj in c(ii:numFeatures)) {
      if(ii!=jj) {
        # a symetric matrix ...
        pairwiseMat[ii,jj] = mine(ecoliData[[ii]], ecoliData[[jj]])$MICe
        pairwiseMat[jj,ii] = mine(ecoliData[[ii]], ecoliData[[jj]])$MICe
      }
    }
  }
  # write to file
  write.table(pairwiseMat, file = micePairwiseFilename, append = FALSE, row.names = F, sep = ",");
}
if(!file.exists(dcorPairwiseFilename)) {
  printf("Computing dCor Pairwise Dependencies for eColi Data")
  pairwiseMat = matrix(0, nrow=numFeatures, ncol=numFeatures)
  for(ii in c(1:numFeatures)) {
    for(jj in c(ii:numFeatures)) {
      if(ii!=jj) {
        # a symetric matrix ...
        pairwiseMat[ii,jj] = dcor(ecoliData[[ii]], ecoliData[[jj]])
        pairwiseMat[jj,ii] = dcor(ecoliData[[ii]], ecoliData[[jj]])
      }
    }
  }
  # write to file
  write.table(pairwiseMat, file = dcorPairwiseFilename, append = FALSE, row.names = F, sep = ",");
}
if(!file.exists(rdcPairwiseFilename)) {
  printf("Computing RDC Pairwise Dependencies for eColi Data")
  pairwiseMat = matrix(0, nrow=numFeatures, ncol=numFeatures)
  for(ii in c(1:numFeatures)) {
    for(jj in c(ii:numFeatures)) {
      if(ii!=jj) {
        # a symetric matrix ...
        pairwiseMat[ii,jj] = rdc(ecoliData[[ii]], ecoliData[[jj]])
        pairwiseMat[jj,ii] = rdc(ecoliData[[ii]], ecoliData[[jj]])
      }
    }
  }
  # write to file
  write.table(pairwiseMat, file = rdcPairwiseFilename, append = FALSE, row.names = F, sep = ",");
}
if(!file.exists(ccorPairwiseFilename)) {
  printf("Computing cCor Pairwise Dependencies for eColi Data")
  pairwiseMat = matrix(0, nrow=numFeatures, ncol=numFeatures)
  for(ii in c(1:numFeatures)) {
    for(jj in c(ii:numFeatures)) {
      if(ii!=jj) {
        # a symetric matrix ...
        pairwiseMat[ii,jj] = rcd(ecoliData[[ii]], ecoliData[[jj]])
        pairwiseMat[jj,ii] = rcd(ecoliData[[ii]], ecoliData[[jj]])
      }
    }
  }
  # write to file
  write.table(pairwiseMat, file = ccorPairwiseFilename, append = FALSE, row.names = F, sep = ",");
}

# read in the "true" network
trueNetwork = read.table("ecoli_results/true.txt",header=F)
trueNetwork = as.matrix(trueNetwork)

# read in computed metrics for the 5 metrics
cosResults = read.csv(cosPairwiseFilename,header=T)
cosResults = as.matrix(cosResults)
miceResults = read.csv(micePairwiseFilename,header=T)
miceResults = as.matrix(miceResults)
dcorResults = read.csv(dcorPairwiseFilename,header=T)
dcorResults = as.matrix(dcorResults)
rdcResults = read.csv(rdcPairwiseFilename,header=T)
rdcResults = as.matrix(rdcResults)
ccorResults = read.csv(ccorPairwiseFilename,header=T)
ccorResults = as.matrix(ccorResults)

# compute ROC and other metrics for each 
mrnet_CosNetwork = mrnet(cosResults)
mrnet_Cos.tbl = validate(mrnet_CosNetwork,trueNetwork)
mrnet_MICeNetwork = mrnet(miceResults)
mrnet_MICe.tbl = validate(mrnet_MICeNetwork,trueNetwork)
mrnet_dCorNetwork = mrnet(dcorResults)
mrnet_dCor.tbl = validate(mrnet_dCorNetwork,trueNetwork)
mrnet_RDCNetwork = mrnet(rdcResults)
mrnet_RDC.tbl = validate(mrnet_RDCNetwork,trueNetwork)
mrnet_cCorNetwork = mrnet(ccorResults)
mrnet_cCor.tbl = validate(mrnet_cCorNetwork,trueNetwork)

# compute ROC
cos_roc = auc.roc(mrnet_Cos.tbl)
mice_roc = auc.roc(mrnet_MICe.tbl)
dcor_roc = auc.roc(mrnet_dCor.tbl)
rdc_roc = auc.roc(mrnet_RDC.tbl)
ccor_roc = auc.roc(mrnet_cCor.tbl)

# compute max(f-score)
cos_maxf = max(fscores(mrnet_Cos.tbl))
mice_maxf = max(fscores(mrnet_MICe.tbl))
dcor_maxf = max(fscores(mrnet_dCor.tbl))
rdc_maxf = max(fscores(mrnet_RDC.tbl))
ccor_maxf = max(fscores(mrnet_cCor.tbl))

# print results
printf("ROC --> CoS=%0.02f MICe=%0.02f dCor=%0.02f RDC=%0.02f cCor=%0.02f",
       cos_roc, mice_roc, dcor_roc, rdc_roc, ccor_roc)
printf("max(F-Score) --> CoS=%0.02f MICe=%0.02f dCor=%0.02f RDC=%0.02f cCor=%0.02f",
       cos_maxf, mice_maxf, dcor_maxf, rdc_maxf, ccor_maxf)

# plot the results
plot(sin, -pi, 2*pi)
dev <- show.roc(mrnet_Cos.tbl,col="dodgerblue",type="s",lwd=3,lty=1,pch=1)
dev <- show.roc(mrnet_dCor.tbl,col="darkorange",device=dev,type="s",lwd=3,lty=2,pch=2)
dev <- show.roc(mrnet_MICe.tbl,col="chartreuse",device=dev,type="s",lwd=3,lty=3,pch=3)

dev2 <- show.roc(mrnet_RDC.tbl,col="firebrick1", type="s",lwd=3,lty=4,pch=4)
dev2 <- show.roc(mrnet_cCor.tbl,col="mediumpurple",device=dev2,type="s",lwd=3,lty=5,pch=5)
#legend("bottomright",c("CoS","dCor","MICe", "RDC", "cCor"), pch = c('-','-','-','-','-'), 
#       col = c("dodgerblue","darkorange","chartreuse", "firebrick1", "mediumpurple"))