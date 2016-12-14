# Ripley Forms Testing and plotting

cat("\014")
source("algorithms.R")
source("copulastatistic.R")
source("multiplot.R")
library(VineCopula)

library(ggplot2)

cat("\014")

set.seed(123)

M = 1000;
form1 <- lcg(65, 1, 2048, M);
form2 <- lcg(1229, 1, 2048, M);
form3 <- lcg(5, 1, 2048, M);
form4 <- lcg(129, 1, 2^64, M);

# plot the returns against each other
p1 <- ggplot(form1, aes(x=x, y=y)) + geom_point(alpha=0.4, colour = "blue", size = 2) + 
  theme(panel.background = element_blank(), 
        axis.title.x = element_text(face="bold", size=13), 
        axis.title.y = element_text(face="bold", size=13), 
        axis.text.x = element_text(size=11), 
        axis.text.y = element_text(size=11), 
        panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  labs(title="Form 1")
p3 <- ggplot(form2, aes(x=x, y=y))   + geom_point(alpha=0.4, colour = "blue", size = 2) + 
  theme(panel.background = element_blank(), 
        axis.title.x = element_text(face="bold", size=13), 
        axis.title.y = element_text(face="bold", size=13), 
        axis.text.x = element_text(size=11), 
        axis.text.y = element_text(size=11),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  labs(title="Form 2")
p2 <- ggplot(form3, aes(x=x, y=y)) + geom_point(alpha=0.4, colour = "blue", size = 2) + 
  theme(panel.background = element_blank(), 
        axis.title.x = element_text(face="bold", size=13), 
        axis.title.y = element_text(face="bold", size=13), 
        axis.text.x = element_text(size=11), 
        axis.text.y = element_text(size=11),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +  
  labs(title="Form 3")
p4 <- ggplot(form4, aes(x=x, y=y)) + geom_point(alpha=0.4, colour = "blue", size = 2) + 
  theme(panel.background = element_blank(), 
        axis.title.x = element_text(face="bold", size=13),
        axis.title.y = element_text(face="bold", size=13),
        axis.text.x = element_text(size=11), 
        axis.text.y = element_text(size=11),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(title="Form 4")

multiplot(p1, p2, p3, p4, cols=2)

# perform monte-carlo simulations to get the average value of each metric for
# these kinds of dependencies

nsim = 1000;
M = 1000;

form1_res = c(0,0,0,0,0)
form2_res = c(0,0,0,0,0)
form3_res = c(0,0,0,0,0)
form4_res = c(0,0,0,0,0)
copu1_res = c(0,0,0,0,0)
copu2_res = c(0,0,0,0,0)
copu3_res = c(0,0,0,0,0)
copu4_res = c(0,0,0,0,0)
copu5_res = c(0,0,0,0,0)

for(ii in c(1:nsim)) {
  printf("sim=%d/1000",ii)
  
  form1 <- lcg(65, 1, 2048, M);
  form2 <- lcg(1229, 1, 2048, M);
  form3 <- lcg(5, 1, 2048, M);
  form4 <- lcg(129, 1, 2^64, M);
  # generate gaussian copula data
  cop1dat <- rCopula(M, normalCopula(0.1))
  # generate gumbel copula data
  cop2dat <- rCopula(M, gumbelCopula(5))
  # generate clayton copula data
  cop3dat <- rCopula(M, claytonCopula(-0.88))
  # generate galambos copula data
  cop4dat <- rCopula(M, galambosCopula(2))
  # generate bb6 copula data
  cop5dat <- rCopula(M, BB6Copula(param = c(2, 2)))
  
  # compute the values for cos,dcor,mice,rdc,cCor for all of these ripley forms
  form1_res[1] = form1_res[1] + cosdv(form1$x, form1$y)
  form1_res[2] = form1_res[2] + dcor(form1$x, form1$y)
  form1_res[3] = form1_res[3] + mine(form1$x, form1$y)$MICe
  form1_res[4] = form1_res[4] + rdc(form1$x, form1$y)
  form1_res[5] = form1_res[5] + rcd(form1$x, form1$y)
   
  
  form2_res[1] = form2_res[1] + cosdv(form2$x, form2$y)
  form2_res[2] = form2_res[2] + dcor(form2$x, form2$y)
  form2_res[3] = form2_res[3] + mine(form2$x, form2$y)$MICe
  form2_res[4] = form2_res[4] + rdc(form2$x, form2$y)
  form2_res[5] = form2_res[5] + rcd(form2$x, form2$y)
  
  form3_res[1] = form3_res[1] + cosdv(form3$x, form3$y)
  form3_res[2] = form3_res[2] + dcor(form3$x, form3$y)
  form3_res[3] = form3_res[3] + mine(form3$x, form3$y)$MICe
  form3_res[4] = form3_res[4] + rdc(form3$x, form3$y)
  form3_res[5] = form3_res[5] + rcd(form3$x, form3$y)
  
  form4_res[1] = form4_res[1] + cosdv(form4$x, form4$y)
  form4_res[2] = form4_res[2] + dcor(form4$x, form4$y)
  form4_res[3] = form4_res[3] + mine(form4$x, form4$y)$MICe
  form4_res[4] = form4_res[4] + rdc(form4$x, form4$y)
  form4_res[5] = form4_res[5] + rcd(form4$x, form4$y)
  
  copu1_res[1] = copu1_res[1] + cosdv(cop1dat[,1], cop1dat[,2])
  copu1_res[2] = copu1_res[2] + dcor(cop1dat[,1], cop1dat[,2])
  copu1_res[3] = copu1_res[3] + mine(cop1dat[,1], cop1dat[,2])$MICe
  copu1_res[4] = copu1_res[4] + rdc(cop1dat[,1], cop1dat[,2])
  copu1_res[5] = copu1_res[5] + rcd(cop1dat[,1], cop1dat[,2])
  
  copu2_res[1] = copu2_res[1] + cosdv(cop2dat[,1], cop2dat[,2])
  copu2_res[2] = copu2_res[2] + dcor(cop2dat[,1], cop2dat[,2])
  copu2_res[3] = copu2_res[3] + mine(cop2dat[,1], cop2dat[,2])$MICe
  copu2_res[4] = copu2_res[4] + rdc(cop2dat[,1], cop2dat[,2])
  copu2_res[5] = copu2_res[5] + rcd(cop2dat[,1], cop2dat[,2])
  
  copu3_res[1] = copu3_res[1] + cosdv(cop3dat[,1], cop3dat[,2])
  copu3_res[2] = copu3_res[2] + dcor(cop3dat[,1], cop3dat[,2])
  copu3_res[3] = copu3_res[3] + mine(cop3dat[,1], cop3dat[,2])$MICe
  copu3_res[4] = copu3_res[4] + rdc(cop3dat[,1], cop3dat[,2])
  copu3_res[5] = copu3_res[5] + rcd(cop3dat[,1], cop3dat[,2])
  
  copu4_res[1] = copu4_res[1] + cosdv(cop4dat[,1], cop4dat[,2])
  copu4_res[2] = copu4_res[2] + dcor(cop4dat[,1], cop4dat[,2])
  copu4_res[3] = copu4_res[3] + mine(cop4dat[,1], cop4dat[,2])$MICe
  copu4_res[4] = copu4_res[4] + rdc(cop4dat[,1], cop4dat[,2])
  copu4_res[5] = copu4_res[5] + rcd(cop4dat[,1], cop4dat[,2])
  
  copu5_res[1] = copu5_res[1] + cosdv(cop5dat[,1], cop5dat[,2])
  copu5_res[2] = copu5_res[2] + dcor(cop5dat[,1], cop5dat[,2])
  copu5_res[3] = copu5_res[3] + mine(cop5dat[,1], cop5dat[,2])$MICe
  copu5_res[4] = copu5_res[4] + rdc(cop5dat[,1], cop5dat[,2])
  copu5_res[5] = copu5_res[5] + rcd(cop5dat[,1], cop5dat[,2])
}

# compute the averages
form1_res = form1_res/nsim
form2_res = form2_res/nsim
form3_res = form3_res/nsim
form4_res = form4_res/nsim
copu1_res = copu1_res/nsim
copu2_res = copu2_res/nsim
copu3_res = copu3_res/nsim
copu4_res = copu4_res/nsim
copu5_res = copu5_res/nsim

printf("Form 1 --> CoS=%0.02f dCor=%0.02f MICe=%0.02f RDC=%0.02f cCor=%0.02f", 
       form1_res[1], form1_res[2], form1_res[3], form1_res[4], form1_res[5])
printf("Form 2 --> CoS=%0.02f dCor=%0.02f MICe=%0.02f RDC=%0.02f cCor=%0.02f", 
       form2_res[1], form2_res[2], form2_res[3], form2_res[4], form2_res[5])
printf("Form 3 --> CoS=%0.02f dCor=%0.02f MICe=%0.02f RDC=%0.02f cCor=%0.02f", 
       form3_res[1], form3_res[2], form3_res[3], form3_res[4], form3_res[5])
printf("Form 4 --> CoS=%0.02f dCor=%0.02f MICe=%0.02f RDC=%0.02f cCor=%0.02f", 
       form4_res[1], form4_res[2], form4_res[3], form4_res[4], form4_res[5])
printf("Gaussian Copula --> CoS=%0.02f dCor=%0.02f MICe=%0.02f RDC=%0.02f cCor=%0.02f", 
       copu1_res[1], copu1_res[2], copu1_res[3], copu1_res[4], copu1_res[5])
printf("Gumbel Copula --> CoS=%0.02f dCor=%0.02f MICe=%0.02f RDC=%0.02f cCor=%0.02f", 
       copu2_res[1], copu2_res[2], copu2_res[3], copu2_res[4], copu2_res[5])
printf("Clayton Copula --> CoS=%0.02f dCor=%0.02f MICe=%0.02f RDC=%0.02f cCor=%0.02f", 
       copu3_res[1], copu3_res[2], copu3_res[3], copu3_res[4], copu3_res[5])
printf("Galambos Copula --> CoS=%0.02f dCor=%0.02f MICe=%0.02f RDC=%0.02f cCor=%0.02f", 
       copu4_res[1], copu4_res[2], copu4_res[3], copu4_res[4], copu4_res[5])
printf("BB6 Copula --> CoS=%0.02f dCor=%0.02f MICe=%0.02f RDC=%0.02f cCor=%0.02f", 
       copu5_res[1], copu5_res[2], copu5_res[3], copu5_res[4], copu5_res[5])
