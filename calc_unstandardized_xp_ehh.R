library(plyr)  ## requires plyr package that contains count() function

################
## calc_unstandardized_xp_ehh.R
##
## @authors Martyna Lukaszewicz
## @contact martyna@uidaho.edu
################
## @description Simulate unstandardized XP-EHH. 
##              Formula described in A Map of Recent Positive Selection in the Human Genome
##              by Voight et al., eq (1). Webpage to the article:
##              https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0040072 
##              
## @param popX : population X in a data frame format, with rows corresponding to individuals in a population,
##               and columns corresponding to bi-allelelic SNPs.
## @param popX : population X in a data frame format, with rows corresponding to individuals in a population,
##               and columns corresponding to bi-allelelic SNPs.
##
## @lastChange 2019-03-01
##
## @changes
##  
##
## Example usage:
##
## Input 2 populations, popX and PopY, in a data frame format, with rows corresponding to individuals in a population,
## and columns corresponding to bi-allelelic SNPs. In the example each population is of size 4, with 6 SNPs per individual.
##
## popX <- as.data.frame(matrix(c(1,1,0,0,0,1,1,0,1,0,1,1,0,1,0,0,0,0,1,1,0,1,0,1),nrow=4,ncol=6,byrow=T))
## popY <-as.data.frame(matrix(c(0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,1,0,0,0,0,0),nrow=4,ncol=6,byrow=T))
## 
###############################
calc_unstandardized_xp_ehh <- function(popX,popY){
  n.col <- ncol(popX)  ## assumed ncol(popX)=ncol(popY)
  n.row <- nrow(popX)  ## assumed nrow(popX)=nrow(popY)
  colnames(popX) <- 1:n.col
  colnames(popY) <- 1:n.col
  ## can ignore all possible combinations if some of them are observed 0 times because choose(0,2)=0 
  numerator_popX <- matrix(rep(NA,n.col*n.col),nrow=n.col,ncol=n.col)
  denominator_popX <- matrix(rep(NA,n.col*n.col),nrow=n.col,ncol=n.col)
  ehh_popX <- matrix(rep(NA,n.col*n.col),nrow=n.col,ncol=n.col)
  
  numerator_popY <- matrix(rep(NA,n.col*n.col),nrow=n.col,ncol=n.col)
  denominator_popY <- matrix(rep(NA,n.col*n.col),nrow=n.col,ncol=n.col)
  ehh_popY <- matrix(rep(NA,n.col*n.col),nrow=n.col,ncol=n.col)
  
  for (i in 1:n.col){
    J <- i:n.col
    for (j in J){
      ## do EHH calculations for popX
      numerator_popX[i,j] <- sum(choose(count(popX, vars=i:j)[,ncol(count(popX, vars=i:j))],rep(2,nrow(count(popX, vars=i:j)))))  ## complete row filled for SNP #1 position
      n_0_popX <- sum((count(popX, vars=i:j)[,1]==0)*count(popX, vars=i:j)[,ncol(count(popX, vars=i:j))]) 
      n_1_popX <- sum((count(popX, vars=i:j)[,1]==1)*count(popX, vars=i:j)[,ncol(count(popX, vars=i:j))]) 
      denominator_popX[i,j] <- sum(choose(c(n_0_popX,n_1_popX), c(2,2)))
      ehh_popX[i,j] <- numerator_popX[i,j]/denominator_popX[i,j]
      
      ## do EHH calculations for popY
      numerator_popY[i,j] <- sum(choose(count(popY, vars=i:j)[,ncol(count(popY, vars=i:j))],rep(2,nrow(count(popY, vars=i:j)))))  ## complete row filled for SNP #1 position
      n_0_popY <- sum((count(popY, vars=i:j)[,1]==0)*count(popY, vars=i:j)[,ncol(count(popY, vars=i:j))]) 
      n_1_popY <- sum((count(popY, vars=i:j)[,1]==1)*count(popY, vars=i:j)[,ncol(count(popY, vars=i:j))]) 
      denominator_popY[i,j] <- sum(choose(c(n_0_popY,n_1_popY), c(2,2)))
      ehh_popY[i,j] <- numerator_popY[i,j]/denominator_popY[i,j]
    }
    if (i > 1){  ## fill out the lower triangle of the numerator matrix
      K <- 1:(i-1)
      for (k in K){
        ## do EHH calculations for popX
        numerator_popX[i,k] <- sum(choose(count(popX, vars=k:i)[,ncol(count(popX, vars=k:i))],rep(2,nrow(count(popX, vars=k:i))))) 
        n_0_popX <- sum((count(popX, vars=k:i)[,1]==0)*count(popX, vars=k:i)[,ncol(count(popX, vars=k:i))]) 
        n_1_popX <- sum((count(popX, vars=k:i)[,1]==1)*count(popX, vars=k:i)[,ncol(count(popX, vars=k:i))]) 
        denominator_popX[i,k] <- sum(choose(c(n_0_popX,n_1_popX), c(2,2)))
        ehh_popX[i,k] <- numerator_popX[i,k]/denominator_popX[i,k]
        
        ## do EHH calculations for popY
        numerator_popY[i,k] <- sum(choose(count(popY, vars=k:i)[,ncol(count(popY, vars=k:i))],rep(2,nrow(count(popY, vars=k:i))))) 
        n_0_popY <- sum((count(popY, vars=k:i)[,1]==0)*count(popY, vars=k:i)[,ncol(count(popY, vars=k:i))]) 
        n_1_popY <- sum((count(popY, vars=k:i)[,1]==1)*count(popY, vars=k:i)[,ncol(count(popY, vars=k:i))])
        denominator_popY[i,k] <- sum(choose(c(n_0_popY,n_1_popY), c(2,2)))
        ehh_popY[i,k] <- numerator_popY[i,k]/denominator_popY[i,k]
      }
    }
  }
  
  ## Calculate unstandardized XP-EHH between popX and popY (A Map of Recent Positive Selection in the Human Genome
  ##by Voight et al., equation (1)):
  unstandardized_xp_ehh <- log(rowSums(ehh_popX)/rowSums(ehh_popY)) ## XP-EHH vector of the length equal to the number of SNPs,
  ## can ignore x-distance units as they cancel out of the equation.
  ## unstandardized_xp_ehh is approx. 0 when the rate of EHH decay is similar on popX and popY alleles.
  return(unstandardized_xp_ehh)
}
#######################

calc_unstandardized_xp_ehh(popX,popY)
