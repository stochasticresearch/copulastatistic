###########
# An R Script which reads in monthly stock prices from 1991-2016 for:
#  SP500
#  Nikkei225
#  DAX
#  FTSE
# It then does a single differencing to calculate "returns" from price,
# performs a Dickey-Fuller test of stationarity, and then applies the different
# dependency metrics on the returns data
###########

cat("\014")
source("algorithms.R")
source("copulastatistic.R")
source("multiplot.R")

library(tseries)
library(ggplot2)
library(gridExtra)

# Read in the data :)
sp500.data <- read.csv(file="./financial_data/SP500.csv",head=TRUE,sep=",")
dax.data   <- read.csv(file="./financial_data/DAX.csv",head=TRUE,sep=",")
nk225.data <- read.csv(file="./financial_data/N225.csv",head=TRUE,sep=",")
cac40.data <- read.csv(file="./financial_data/CAC.csv",head=TRUE, sep=",")

# we only care about closing price
sp500.close = sp500.data$Close
dax.close = dax.data$Close
nk225.close = nk225.data$Close
cac40.close = cac40.data$Close

sp500.returns = sp500.close[2:length(sp500.close)]-sp500.close[1:length(sp500.close)-1]
dax.returns   = dax.close[2:length(dax.close)]-dax.close[1:length(dax.close)-1]
nk225.returns = nk225.close[2:length(nk225.close)]-nk225.close[1:length(nk225.close)-1]
cac40.returns = cac40.close[2:length(cac40.close)]-cac40.close[1:length(cac40.close)-1]

returns.data = data.frame("SP500"=sp500.returns,"DAX"=dax.returns,"NK225"=nk225.returns,"CAC40"=cac40.returns)

# apply Dickey-Fuller test (if p<alpha, then we reject H0 that a unit-root exist, 
#                           meaning the sequence is stationary)
sp500_dfTest <- adf.test(sp500.returns)
dax_dfTest <- adf.test(dax.returns)
nk225_dfTest <- adf.test(nk225.returns)
cac40_dfTest <- adf.test(cac40.returns)

# apply Philips-Perron test (if p<alpha, then we reject H0 that a unit-root exist, 
#                            meaning the sequence is stationary)
sp500_ppTest <- pp.test(sp500.returns)
dax_ppTest   <- pp.test(dax.returns)
nk225_ppTest <- pp.test(nk225.returns)
cac40_ppTest <- pp.test(cac40.returns)

# plot the returns against each other
p1 <- ggplot(returns.data, aes(x=SP500, y=NK225)) + geom_point(alpha=0.4, colour = "blue", size = 2) + 
      theme(panel.background = element_blank(), 
            panel.grid.major = element_line(colour="black", size=0.15), 
            panel.grid.minor = element_blank(), 
            axis.title.x = element_text(face="bold", size=13), 
            axis.title.y = element_text(face="bold", size=13), 
            axis.text.x = element_text(size=11), 
            axis.text.y = element_text(size=11), 
            panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
      scale_x_continuous(breaks = seq(-200, 200, 100)) + scale_y_continuous(breaks = seq(-2500, 2500, 1500)) + 
      labs(x="SP500 Returns", y="NK225 Returns")
p2 <- ggplot(returns.data, aes(x=CAC40, y=DAX))   + geom_point(alpha=0.4, colour = "blue", size = 2) + 
      theme(panel.background = element_blank(), 
            panel.grid.major = element_line(colour="black", size=0.15), 
            panel.grid.minor = element_blank(), 
            axis.title.x = element_text(face="bold", size=13), 
            axis.title.y = element_text(face="bold", size=13), 
            axis.text.x = element_text(size=11), 
            axis.text.y = element_text(size=11),
            panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
      scale_x_continuous(breaks = seq(-500, 900, 400)) + scale_y_continuous(breaks = seq(-1350, 1600, 750)) + 
      labs(x="CAC40 Returns", y="DAX Returns")
p3 <- ggplot(returns.data, aes(x=SP500, y=CAC40)) + geom_point(alpha=0.4, colour = "blue", size = 2) + 
      theme(panel.background = element_blank(), 
            panel.grid.major = element_line(colour="black", size=0.15), 
            panel.grid.minor = element_blank(), 
            axis.title.x = element_text(face="bold", size=13), 
            axis.title.y = element_text(face="bold", size=13), 
            axis.text.x = element_text(size=11), 
            axis.text.y = element_text(size=11),
            panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
      scale_x_continuous(breaks = seq(-125, 200, 100)) + scale_y_continuous(breaks = seq(-400, 800, 500)) + 
      labs(x="SP500 Returns", y="CAC40 Returns")
p4 <- ggplot(returns.data, aes(x=DAX, y=NK225)) + geom_point(alpha=0.4, colour = "blue", size = 2) + 
      theme(panel.background = element_blank(), 
            panel.grid.major = element_line(colour="black", size=0.15), 
            panel.grid.minor = element_blank(), 
            axis.title.x = element_text(face="bold", size=13),
            axis.title.y = element_text(face="bold", size=13),
            axis.text.x = element_text(size=11), 
            axis.text.y = element_text(size=11),
            panel.border = element_rect(colour = "black", fill=NA, size=1)) +
      scale_x_continuous(breaks = seq(-1100, 1600, 750)) + scale_y_continuous(breaks = seq(-2500, 2500, 1500)) + labs(x="DAX Returns", y="NK225 Returns")

multiplot(p1, p2, p3, p4, cols=2)

# compute pair-wise dependencies w/ MIC, CoS, RDC, TDC, dCorr
sp_cac.cos  = cosdv(sp500.returns, cac40.returns)
sp_cac.mic  = mine(sp500.returns, cac40.returns)$MICe
sp_cac.rdc  = rdc(sp500.returns, cac40.returns)
sp_cac.ccor = rcd(sp500.returns, cac40.returns)
sp_cac.dcor = dcor(sp500.returns, cac40.returns)

sp_dax.cos  = cosdv(sp500.returns, dax.returns)
sp_dax.mic  = mine(sp500.returns, dax.returns)$MICe
sp_dax.rdc  = rdc(sp500.returns, dax.returns)
sp_dax.ccor = rcd(sp500.returns, dax.returns)
sp_dax.dcor = dcor(sp500.returns, dax.returns)

sp_nk.cos  = cosdv(sp500.returns, nk225.returns)
sp_nk.mic  = mine(sp500.returns, nk225.returns)$MICe
sp_nk.rdc  = rdc(sp500.returns, nk225.returns)
sp_nk.ccor = rcd(sp500.returns, nk225.returns)
sp_nk.dcor = dcor(sp500.returns, nk225.returns)

cac_dax.cos  = cosdv(cac40.returns, dax.returns)
cac_dax.mic  = mine(cac40.returns, dax.returns)$MICe
cac_dax.rdc  = rdc(cac40.returns, dax.returns)
cac_dax.ccor = rcd(cac40.returns, dax.returns)
cac_dax.dcor = dcor(cac40.returns, dax.returns)

cac_nk.cos  = cosdv(cac40.returns, nk225.returns)
cac_nk.mic  = mine(cac40.returns, nk225.returns)$MICe
cac_nk.rdc  = rdc(cac40.returns, nk225.returns)
cac_nk.ccor = rcd(cac40.returns, nk225.returns)
cac_nk.dcor = dcor(cac40.returns, nk225.returns)

dax_nk.cos  = cosdv(dax.returns, nk225.returns)
dax_nk.mic  = mine(dax.returns, nk225.returns)$MICe
dax_nk.rdc  = rdc(dax.returns, nk225.returns)
dax_nk.ccor = rcd(dax.returns, nk225.returns)
dax_nk.dcor = dcor(dax.returns, nk225.returns)

# compute 3-D financial data analysis
sp_nk_dax = cosdv3d(sp500.returns, nk225.returns, dax.returns)
sp_cac_dax = cosdv3d(sp500.returns, cac40.returns, dax.returns)
sp_cac_nk  = cosdv3d(sp500.returns, cac40.returns, nk225.returns)
cac_nk_dax = cosdv3d(cac40.returns, nk225.returns, dax.returns)

# compute 4-D financial data analysis
sp_cac_nk_dax = cosdv4d(sp500.returns, cac40.returns, nk225.returns, dax.returns)