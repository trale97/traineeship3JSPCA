#source("C:/Users/TSB-MTO/Documents/Tra/JointSPCA_finalizing.R")
#source("C:/Users/TSB-MTO/Documents/Tra/IS_JSPCA.R")

#2 groups and 2 components
#primary loading shifts
#no of loading differences = 4
source("C:/Users/TSB-MTO/Documents/Tra/JointSPCA_finalizing.R")
source("C:/Users/TSB-MTO/Documents/Tra/IS_JSPCA.R")
library(MASS)
n <- c(50, 200, 600, 1000)
final_result <- data.frame()
set.seed(45468)
for (i in 1:length(n)){
  getValue <- function(){
    p1g1 <- c(0,rep(sqrt(.6),9),rep(0,10))
    p2g1 <- c(sqrt(.6),rep(0,9),rep(sqrt(.6),10))
    P1 <- cbind(p1g1,p2g1)
    p1g2 <- c(rep(sqrt(.6),2),0,rep(sqrt(.6),7),rep(0,10))
    p2g2 <- c(0,0,sqrt(.6),rep(0,7),rep(sqrt(.6),10))
    P2 <- cbind(p1g2,p2g2)
    PSI1 <- diag(1-rowSums(P1^2))
    SIGMA1 <- P1%*%t(P1)+PSI1 
    data1 <- mvrnorm(n = n[i], mu = rep(0,20), Sigma = SIGMA1, empirical = TRUE)
    PSI2 <- diag(1-rowSums(P2^2))
    SIGMA2 <- P2%*%t(P2)+PSI2
    data2 <- mvrnorm(n = n[i], mu = rep(0,20), Sigma = SIGMA2, empirical = TRUE)
    DATA <- list(data1,data2)

    #RES1 <-
      #JSPCA(
       # DATA,
       # R=2,
        #TYPE = 'sca',
        #lambda = 0.1,
      #  CARD = c(9,11),
       # MaxIter = 20,
       # eps = 10^-6 
      #)
    
    #RES1$Lossvec
    #plot(RES1$Lossvec[2:length(RES1$Lossvec)])
    #RES1$iterinner
    #RES1$iterouter
    #print the loadings on the 2 components of each group
    #cbind(round(P1,3),round(RES1$loadings[[1]],3),round(P2,3),round(RES1$loadings[[2]],3))
    #p1 is true loading for first group, p2 for 2nd group
    
    lambdavec <- c(0.01,0.05,0.1,0.2,0.5)#index of sparseness for all combinations
    cardvec <- c(10,20,30)#cardinality constraint over all components
    #cardmat <- matrix(rep(c(5,10,15),2),nrow = length(cardvec))#per component
    avec <- matrix(nrow=length(lambdavec)*length(cardvec),ncol=6)
    colnames(avec) <- c('lambda','K',"IS","PEV","Prop0","PropUnique")
    INIT <- 'sca'
    eps <- 10^-6
    
    teller <- 1
    for (j in 1:length(lambdavec)){
      for (k in 1:length(cardvec)){
        a <- IS_JSPCA(DATA, R=2, INIT, lambda = lambdavec[j], card = cardvec[k], MaxIter = 20, eps)
        avec[teller,] <- c(lambdavec[j],cardvec[k],a$value,a$vaf,a$propzero,a$propunique)
        teller <- teller+1
      }
    }
    
    round(avec,3)
    selmodelindex <- which.max(avec[,3]) #best solution
    selmodelindex
    selmodel <- jointSPCA(DATA, R=2, INIT, lambda = avec[selmodelindex,1], CARD = avec[selmodelindex,2], MaxIter = 200, eps)
    cbind(round(P1,3),round(selmodel$loadings[[1]],3),round(P2,3),round(selmodel$loadings[[2]],3))
    
    #calculation zero/non-zero recovery and tucker congruence component scores 
    perm <- gtools::permutations(2, 2)
    corrate1 <- vector(length = nrow(perm))
    corrate2 <- vector(length = nrow(perm))
    spcr_tucongrT1 <- vector(length = nrow(perm))
    spcr_tucongrT2 <- vector(length = nrow(perm))
    for (p in 1:nrow(perm)) {
      corrate1[p]<-num_correct(P1, round(selmodel$loadings[[1]][,perm[p,]]))
      corrate2[p]<-num_correct(P2, round(selmodel$loadings[[2]][,perm[p,]]))
      spcr_tucongrT1[p] <- sum(diag(abs(psych::factor.congruence(P1,selmodel$loadings[[1]][,perm[p,]]))))/2
      spcr_tucongrT2[p] <- sum(diag(abs(psych::factor.congruence(P2,selmodel$loadings[[2]][,perm[p,]]))))/2
      }
    corrate <- mean(max(corrate1),max(corrate2))
    spcr_tucongrT <- mean(max(spcr_tucongrT1),max(spcr_tucongrT2))
    
    parameters <- c(corrate,spcr_tucongrT)
    return(parameters)
  }
  
  result <- as.data.frame(t(replicate(50,getValue())))
  result$size <- c(rep(n[i],50))
  final_result <- rbind(final_result,result)
}

colnames(final_result) <- c("Recovery","Tucker","Group_size")
save(final_result, file = "final_result1.Rda")
