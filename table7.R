# Regenerate table-7

cat("\014")
source("algorithms.R")
source("copulastatistic.R")

set.seed(123)

n = 1000;
nsim = 1000;

p = c(0.5, 1, 2, 3, 4, 5)

for(pval in p) {
  printf("*****************************")
  printf("p = %0.02f", pval)
  
  y1_rdc_total = 0;
  y1_cCor_total = 0;
  y1_cos_total = 0;
  y1_dcor_total = 0;
  y1_mice_total = 0;
  
  y2_rdc_total = 0;
  y2_cCor_total = 0;
  y2_cos_total = 0;
  y2_dcor_total = 0;
  y2_mice_total = 0;
  
  y3_rdc_total = 0;
  y3_cCor_total = 0;
  y3_cos_total = 0;
  y3_dcor_total = 0;
  y3_mice_total = 0;
  
  y4_rdc_total = 0;
  y4_cCor_total = 0;
  y4_cos_total = 0;
  y4_dcor_total = 0;
  y4_mice_total = 0;
  
  for(ii in 1:nsim) {
    x = runif(n, -5, 5);
    
    pe = pval*rnorm(n);
    
    y1 = 2*x + 1; y1 = y1*(1+pe);
    y2 = (x^2-0.25)*(x^2-1); y2 = y2*(1+pe);
    y3 = cos(x); y3 = y3*(1+pe);
    x2 = runif(n, 0, 1);
    y4 = (2*rbinom(n,1,0.5)-1) * (sqrt(1 - (2*x2 - 1)^2)); y4 = y4*(1+pe);
    
    rdcVal = rdc(x,y1);
    cCorVal = rcd(x,y1);
    cosVal = cosdv(x,y1);
    miceVal = mine(x,y1)$MICe;
    dCorVal = dcor(x,y1);
    y1_rdc_total = y1_rdc_total + rdcVal;
    y1_cCor_total = y1_cCor_total + cCorVal;
    y1_cos_total = y1_cos_total + cosVal;
    y1_mice_total = y1_mice_total + miceVal;
    y1_dcor_total = y1_dcor_total + dCorVal;
    
    rdcVal = rdc(x,y2);
    cCorVal = rcd(x,y2);
    cosVal = cosdv(x,y2);
    miceVal = mine(x,y2)$MICe;
    dCorVal = dcor(x,y2);
    y2_rdc_total = y2_rdc_total + rdcVal;
    y2_cCor_total = y2_cCor_total + cCorVal;
    y2_cos_total = y2_cos_total + cosVal;
    y2_mice_total = y2_mice_total + miceVal;
    y2_dcor_total = y2_dcor_total + dCorVal;
    
    rdcVal = rdc(x,y3);
    cCorVal = rcd(x,y3);
    cosVal = cosdv(x,y3);
    miceVal = mine(x,y3)$MICe;
    dCorVal = dcor(x,y3);
    y3_rdc_total = y3_rdc_total + rdcVal;
    y3_cCor_total = y3_cCor_total + cCorVal;
    y3_cos_total = y3_cos_total + cosVal;
    y3_mice_total = y3_mice_total + miceVal;
    y3_dcor_total = y3_dcor_total + dCorVal;
    
    rdcVal = rdc(x2,y4);
    cCorVal = rcd(x2,y4);
    cosVal = cosdv(x2,y4);
    miceVal = mine(x2,y4)$MICe;
    dCorVal = dcor(x2,y4);
    y4_rdc_total = y4_rdc_total + rdcVal;
    y4_cCor_total = y4_cCor_total + cCorVal;
    y4_cos_total = y4_cos_total + cosVal;
    y4_mice_total = y4_mice_total + miceVal;
    y4_dcor_total = y4_dcor_total + dCorVal;
  }
  
  y1_rdc_av  = y1_rdc_total/nsim
  y1_cCor_av = y1_cCor_total/nsim
  y1_cos_av  = y1_cos_total/nsim
  y1_mice_av = y1_mice_total/nsim
  y1_dcor_av = y1_dcor_total/nsim
  
  y2_rdc_av  = y2_rdc_total/nsim
  y2_cCor_av = y2_cCor_total/nsim
  y2_cos_av  = y2_cos_total/nsim
  y2_mice_av = y2_mice_total/nsim
  y2_dcor_av = y2_dcor_total/nsim
  
  y3_rdc_av  = y3_rdc_total/nsim
  y3_cCor_av = y3_cCor_total/nsim
  y3_cos_av  = y3_cos_total/nsim
  y3_mice_av = y3_mice_total/nsim
  y3_dcor_av = y3_dcor_total/nsim
  
  y4_rdc_av  = y4_rdc_total/nsim
  y4_cCor_av = y4_cCor_total/nsim
  y4_cos_av  = y4_cos_total/nsim
  y4_mice_av = y4_mice_total/nsim
  y4_dcor_av = y4_dcor_total/nsim
  
  printf("y1 --> CoS=%0.02f dCor=%0.02f MICe=%0.02f RDC=%0.02f cCor=%0.02f", 
         y1_cos_av, y1_dcor_av, y1_mice_av, y1_rdc_av, y1_cCor_av)
  printf("y2 --> CoS=%0.02f dCor=%0.02f MICe=%0.02f RDC=%0.02f cCor=%0.02f", 
         y2_cos_av, y2_dcor_av, y2_mice_av, y2_rdc_av, y2_cCor_av)
  printf("y3 --> CoS=%0.02f dCor=%0.02f MICe=%0.02f RDC=%0.02f cCor=%0.02f", 
         y3_cos_av, y3_dcor_av, y3_mice_av, y3_rdc_av, y3_cCor_av)
  printf("y4 --> CoS=%0.02f dCor=%0.02f MICe=%0.02f RDC=%0.02f cCor=%0.02f", 
         y4_cos_av, y4_dcor_av, y4_mice_av, y4_rdc_av, y4_cCor_av)
  printf("*****************************")
}

