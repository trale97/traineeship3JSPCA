######## MGEFA by De Roover and Vermunt (2019) #################
## Author: Tra Le
##################### Preliminary ##############################

## 2 groups 2 components
## PL shift 
## 4 differences
# load library
rm(list = ls())
setwd("C:/Users/onno1234/Documents/MGFR")
library(MASS)
library(haven)
LG <- "C:/Users/onno1234/Documents/LatentGOLD6.0/lg60.exe"
# function to generate data 
p1g1 <- c(0,rep(sqrt(.6),9),rep(0,10))
p2g1 <- c(sqrt(.6),rep(0,9),rep(sqrt(.6),10))
P1 <- cbind(p1g1,p2g1)
p1g2 <- c(rep(sqrt(.6),2),0,rep(sqrt(.6),7),rep(0,10))
p2g2 <- c(0,0,sqrt(.6),rep(0,7),rep(sqrt(.6),10))
P2 <- cbind(p1g2,p2g2)
data_generation <- function(N){
  p1g1 <- c(0,rep(sqrt(.6),9),rep(0,10))
  p2g1 <- c(sqrt(.6),rep(0,9),rep(sqrt(.6),10))
  P1 <- cbind(p1g1,p2g1)
  p1g2 <- c(rep(sqrt(.6),2),0,rep(sqrt(.6),7),rep(0,10))
  p2g2 <- c(0,0,sqrt(.6),rep(0,7),rep(sqrt(.6),10))
  P2 <- cbind(p1g2,p2g2)
  PSI1 <- diag(1-rowSums(P1^2))
  SIGMA1 <- P1%*%t(P1)+PSI1 
  data1 <- mvrnorm(n = N, mu = rep(0,20), Sigma = SIGMA1, empirical = TRUE)
  PSI2 <- diag(1-rowSums(P2^2))
  SIGMA2 <- P2%*%t(P2)+PSI2
  data2 <- mvrnorm(n = N, mu = rep(0,20), Sigma = SIGMA2, empirical = TRUE)
  DATA <- list(data1,data2)
}

n_dataset <- 50 #number of datasets per condition
for (k in 1:n_dataset){
  data_list <- data_generation(N=50)
  data_df <- as.data.frame(rbind(data_list[[1]], data_list[[2]]))
  data_df$G <- c(rep(1,50), rep(2,50))
  write_sav(data_df, paste0("data_N50sim",k,"_G",G,"R",R,"PLshift_4diff.sav"))
  print(k)
}


## write LG syntax
# syntax for MGEFA no rotation
NoRotation <- function(syntaxName, infile, outfile, outfile1){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
//LG6.0//
version = 6.0
infile '", infile,"'

model
options
   algorithm 
      tolerance=1e-008 emtolerance=0.01 emiterations=2500 nriterations=500;
   startvalues
      seed=0 sets=5 tolerance=1e-005 iterations=100 PCA;
   bayes
      categorical=0 variances=0 latent=0 poisson=0;
   montecarlo
      seed=0 replicates=500 tolerance=1e-008;
   quadrature  nodes=10;
   missing  excludeall;
   output      
      iterationdetail classification parameters=effect standarderrors probmeans=posterior profile bivariateresiduals writeparameters='", outfile1,"' write='", outfile,"';
variables
   groupid G;
   dependent V1 continuous, V2 continuous, V3 continuous, V4 continuous, V5 continuous, V6 continuous, V7 continuous, V8 continuous, V9 continuous, V10 continuous, V11 continuous, V12 continuous, V13 continuous, V14 continuous, V15 continuous, V16 continuous, V17 continuous, V18 continuous, V19 continuous, V20 continuous;
   independent G nominal;	
   latent
      F1 continuous, 
      F2 continuous;
equations
   F1 | G;
   F2 | G;
   V1 - V20 <- F1 | G + F2 | G;
   V1 - V20 | G;
end model
")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}

# LG syntax for rotation oblimin alignment
Rotated_OA <- function(syntaxName, infile, criterion, outfile, startval){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
//LG6.0//
version = 6.0
infile '", infile,"' 

model
options
   algorithm 
      tolerance=1e-008 emtolerance=0.01 emiterations=0 nriterations=0;
   startvalues
      seed=0 sets=1 tolerance=1e-005 iterations=0 PCA;
   bayes
      categorical=0 variances=0 latent=0 poisson=0;
   montecarlo
      seed=0 replicates=500 tolerance=1e-008;
   quadrature  nodes=10;
   missing  excludeall;
   rotation oblimin alignment = ", criterion," sets=1;
   output      
      iterationdetail classification parameters=effect standarderrors probmeans=posterior profile bivariateresiduals write='", outfile,"';
variables
   groupid G;
   dependent V1 continuous, V2 continuous, V3 continuous, V4 continuous, V5 continuous, V6 continuous, V7 continuous, V8 continuous, V9 continuous, V10 continuous, V11 continuous, V12 continuous, V13 continuous, V14 continuous, V15 continuous, V16 continuous, V17 continuous, V18 continuous, V19 continuous, V20 continuous;
   independent G nominal;	
   latent
      F1 continuous, 
      F2 continuous;
equations
// factor variances
   F1 | G;
   F2 | G;
   F1 <-> F2 | G;
// regression models for items
   V1 - V20 <- F1 | G + F2 | G;
// item variances 
   V1 - V20 | G;
   '", startval,"'
end model
")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}

# LG syntax for oblimin procrustes
Rotated_OP <- function(syntaxName, infile, criterion, outfile, startval){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
//LG6.0//
version = 6.0
infile '", infile,"' 

model
options
   algorithm 
      tolerance=1e-008 emtolerance=0.01 emiterations=0 nriterations=0;
   startvalues
      seed=0 sets=1 tolerance=1e-005 iterations=0 PCA;
   bayes
      categorical=0 variances=0 latent=0 poisson=0;
   montecarlo
      seed=0 replicates=500 tolerance=1e-008;
   quadrature  nodes=10;
   missing  excludeall;
   rotation oblimin procrustes = ", criterion," sets=1;
   output      
      iterationdetail classification parameters=effect standarderrors probmeans=posterior profile bivariateresiduals write='", outfile,"';
variables
   groupid G;
   dependent V1 continuous, V2 continuous, V3 continuous, V4 continuous, V5 continuous, V6 continuous, V7 continuous, V8 continuous, V9 continuous, V10 continuous, V11 continuous, V12 continuous, V13 continuous, V14 continuous, V15 continuous, V16 continuous, V17 continuous, V18 continuous, V19 continuous, V20 continuous;
   independent G nominal;	
   latent
      F1 continuous, 
      F2 continuous;
equations
// factor variances
   F1 | G;
   F2 | G;
   F1 <-> F2 | G;
// regression models for items
   V1 - V20 <- F1 | G + F2 | G;
// item variances 
   V1 - V20 | G;
   '", startval,"'
end model
")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}

# generate different LG syntax for each rotation criterion
crit <- c(.01, .10, .30, .50, .70)
for (j in crit){
  for (k in 1:n_dataset){
    NoRotation(syntaxName = paste0("NoRotation_N50sim",k,"_G",G,"R",R,"PLshift_4diff"), infile = paste0("data_N50sim",k,"_G",G,"R",R,"PLshift_4diff.sav"), outfile1 = paste0("startingvalues_N50sim",k,"_G",G,"R",R,"PLshift_4diff.txt"), outfile = paste0("results_N50sim",k,"_G",G,"R",R,"PLshift_4diff.csv"))
    Rotated_OA(syntaxName = paste0("Rotated_OA_N50OA",j,"sim",k,"_G",G,"R",R,"PLshift_4diff"), infile = paste0("data_N50sim",k,"_G",G,"R",R,"PLshift_4diff.sav"), outfile = paste0("results_N50OA",j,"sim",k, "_G",G,"R",R,"PLshift_4diff.csv"), startval = paste0("startingvalues_N50sim",k,"_G",G,"R",R,"PLshift_4diff.txt"), criterion = j)
    Rotated_OP(syntaxName = paste0("Rotated_OP_N50OP",j,"sim",k,"_G",G,"R",R,"PLshift_4diff"), infile = paste0("data_N50sim",k,"_G",G,"R",R,"PLshift_4diff.sav"), outfile = paste0("results_N50OP",j,"sim",k, "_G",G,"R",R,"PLshift_4diff.csv"), startval = paste0("startingvalues_N50sim",k,"_G",G,"R",R,"PLshift_4diff.txt"), criterion = j)
    print(c(j,k))
  }
}



# run LG syntax no rotation to get starting values
for (k in 1:n_dataset){
  shell(paste0(LG, " ", "NoRotation_N50sim",k, "_G",G,"R",R,"PLshift_4diff.lgs"," ", "/b"))
  print(k)
}


# run LG syntaxes for different rotation criteria

for (j in crit){
  for (k in 1:n_dataset){
    shell(paste0(LG, " ", "Rotated_OA_N50OA",j,"sim",k,"_G",G,"R",R,"PLshift_4diff.lgs", " ","Rotated_OP_N50OP",j,"sim",k, "_G",G,"R",R,"PLshift_4diff.lgs", " ", "/b"))
    print(c(j,k))
  }
}


#extract the loadings and calculate the recovery rate
parameters <- matrix(NA, nrow = n_dataset, ncol = 2)
result <- data.frame()
final_result <- data.frame()
for (j in crit){
  for (k in 1:n_dataset){
    results <- read.csv(paste0("results_N50OA",j,"sim",k, "_G",G,"R",R,"PLshift_4diff.csv"), header = T, sep = "\t")
    loadings <- strsplit(as.character(results[6,]), ',')
    loadingsF1G1 <- matrix(as.numeric(loadings[[1]][c(21:40)]), ncol = 1)
    loadingsF2G1 <- matrix(as.numeric(loadings[[1]][c(41:60)]), ncol = 1)
    loadingsF1G2 <- matrix(as.numeric(loadings[[1]][c(84:103)]), ncol = 1)
    loadingsF2G2 <- matrix(as.numeric(loadings[[1]][c(104:123)]), ncol = 1)
    loadings_2G <- list(G1 = cbind(loadingsF1G1, loadingsF2G1), G2 = cbind(loadingsF1G2, loadingsF2G2))
    perm <- gtools::permutations(2, 2)
    #corrate1 <- vector(length = nrow(perm))
    #corrate2 <- vector(length = nrow(perm))
    spcr_tucongrT1 <- vector(length = nrow(perm))
    spcr_tucongrT2 <- vector(length = nrow(perm))
    for (p in 1:nrow(perm)) {
      #corrate1[p]<-num_correct(P1, round(selmodel$loadings[[1]][,perm[p,]]))
      #corrate2[p]<-num_correct(P2, round(selmodel$loadings[[2]][,perm[p,]]))
      spcr_tucongrT1[p] <- sum(diag(abs(psych::factor.congruence(P1,loadings_2G$G1[,perm[p,]]))))/2
      spcr_tucongrT2[p] <- sum(diag(abs(psych::factor.congruence(P2,loadings_2G$G2[,perm[p,]]))))/2
    }
    #corrate <- mean(max(corrate1),max(corrate2))
    spcr_tucongrT <- mean(max(spcr_tucongrT1),max(spcr_tucongrT2))
    
    parameters[k,] <- c(spcr_tucongrT,j)
  }
  result <- rbind(result, parameters)
}


R <- 2
G <- 2
save(result, file =paste0("resultOA_G",G,"R",R,"PLshift_4diff.Rda"))

parameters2 <- matrix(NA, nrow = n_dataset, ncol = 2)
result2 <- data.frame()
final_result <- data.frame()
for (j in crit){
  for (k in 1:n_dataset){
    results <- read.csv(paste0("results_N50OP",j,"sim",k, "_G",G,"R",R,"PLshift_4diff.csv"), header = T, sep = "\t")
    loadings <- strsplit(as.character(results[6,]), ',')
    loadingsF1G1 <- matrix(as.numeric(loadings[[1]][c(21:40)]), ncol = 1)
    loadingsF2G1 <- matrix(as.numeric(loadings[[1]][c(41:60)]), ncol = 1)
    loadingsF1G2 <- matrix(as.numeric(loadings[[1]][c(84:103)]), ncol = 1)
    loadingsF2G2 <- matrix(as.numeric(loadings[[1]][c(104:123)]), ncol = 1)
    loadings_2G <- list(G1 = cbind(loadingsF1G1, loadingsF2G1), G2 = cbind(loadingsF1G2, loadingsF2G2))
    perm <- gtools::permutations(2, 2)
    #corrate1 <- vector(length = nrow(perm))
    #corrate2 <- vector(length = nrow(perm))
    spcr_tucongrT1 <- vector(length = nrow(perm))
    spcr_tucongrT2 <- vector(length = nrow(perm))
    for (p in 1:nrow(perm)) {
      #corrate1[p]<-num_correct(P1, round(selmodel$loadings[[1]][,perm[p,]]))
      #corrate2[p]<-num_correct(P2, round(selmodel$loadings[[2]][,perm[p,]]))
      spcr_tucongrT1[p] <- sum(diag(abs(psych::factor.congruence(P1,loadings_2G$G1[,perm[p,]]))))/2
      spcr_tucongrT2[p] <- sum(diag(abs(psych::factor.congruence(P2,loadings_2G$G2[,perm[p,]]))))/2
    }
    #corrate <- mean(max(corrate1),max(corrate2))
    spcr_tucongrT <- mean(max(spcr_tucongrT1),max(spcr_tucongrT2))
    
    parameters2[k,] <- c(spcr_tucongrT,j)
  }
  result2 <- rbind(result2, parameters2)
}
save(result2, file = paste0("resultOP_G",G,"R",R,"PLshift_4diff.Rda"))


