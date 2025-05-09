rm(list=ls())

#load packages
library(ivreg)
library(MASS)
library(glmnet)
library(foreach)
library(doParallel)
library(Matrix)
library(sisVIVE)
library(ncvreg)
library(screening)
library(bestsubset)
library(gurobi)
library(ManyIV)
library(bigstatsr)
library(bigsnpr)

#load relevant functions
source("../all functions.R")

#Input
#G: N * p, matrix, no NA
#X: N * 1, exposure
#Y: N * 1, outcome

set.seed(927)
N=nrow(G)
p=ncol(G)
repeat_time=30#repeat times of MR-SPLIT+
k=2#split the samples into to subsets when applying MR-SPLIT+

#Center
X_c=X-mean(X)
Y_c=Y-mean(Y)
G_c=apply(G, 2, function(col) col - mean(col))
#columns of G scale to unit length
G_unit <- apply(G_c, 2, function(col) col / sqrt(sum(col^2)))
rm(G_c,G)

#SIS  
fit=screening(G_unit,X_c,method = "sis",num.select = 80)
SIS_select=fit$screen
G_candidate=G_unit[,SIS_select]
rm(G_unit)
p_n=ncol(G_candidate)

#MR-SPLIT+
hatX=matrix(0,ncol=repeat_time,nrow=N)#estimated exposure from stage one
MR_relevant=list()#relevant IVs
MR_SPLIT_Results=matrix(0,ncol = 5,nrow = repeat_time)
colnames(MR_SPLIT_Results)=c("Estimation","Std","pvalue","lc","uc")
G_invalid=list()#invalid IVs
#multiple split
for(t in 1:repeat_time){
  #stage one
  MR_stageone=MR_SPLIT_one_CIfirst(X_c,Y_c,G_candidate,k=2,F_threshold=1000000)
  hatX[,t]=MR_stageone$hatX
  
  MR_relevant[[t]]=Reduce(union, MR_stageone$selectIV[1:k])
  
  #Stage two
  # Compute XtX_inv as a scalar
  XtX_inv <- 1 / sum(hatX[, t]^2)  # The reciprocal of the sum of squares of a vector
  
  # calculate the projection part of hatX
  H_vector <- hatX[, t] * XtX_inv  # equivalent to hatX[, t] %*% solve(t(hatX[, t]) %*% hatX[, t])
  
  # initialize the result matrix
  Z_trans <- matrix(0, nrow = nrow(G_candidate), ncol = length(MR_relevant[[t]]))
  
  # process G_candidate column by column
  for (j in seq_along(MR_relevant[[t]])) {
    G_col <- G_candidate[, MR_relevant[[t]][j], drop = FALSE]  # Extract a column
    Z_trans[, j] <- G_col - H_vector * sum(hatX[, t] * G_col)  # projection operation
  }
  
  
  #test if there are any invalid IVs  
  many_est = IVreg(Y_c~-1+X_c|G_candidate[,MR_relevant[[t]]],  inference = "re") ## Depends on ManyIV package
  test = IVoverid(many_est)
  if(test[2,2]<0.05){
    for(ii in 1:(ncol(Z_trans)-1)){
      results=bs(Z_trans,Y_c,k=ii)
      invalid_index=which(results$beta!=0)
      valid_index=which(results$beta==0)
      many_est = IVreg(Y_c~-1+X_c+G_candidate[,MR_relevant[[t]][invalid_index]]|G_candidate[,MR_relevant[[t]][valid_index]],  inference = "re") ## Depends on ManyIV package
      test = IVoverid(many_est)
      if(test[2,2]>0.05)break
    }
    optimal_k=ii
    results=bs(Z_trans,Y_c,k=optimal_k)
    
    G_invalid[[t]]=MR_relevant[[t]][results$beta!=0]
    
  }else{
    G_invalid[[t]]=NA
  }
  
  if(sum(is.na(G_invalid[[t]]))==0){
    fit=lm(Y_c~hatX[,t]+ G_candidate[,G_invalid[[t]]])
  }else{
    fit=lm(Y_c~hatX[,t])
  }
  MR_SPLIT=c(summary(fit)$coef[2,c(1,2,4)])
  MR_SPLIT_Results[t,]=c(MR_SPLIT,MR_SPLIT[1]-qnorm(0.975)*MR_SPLIT[2],MR_SPLIT[1]+qnorm(0.975)*MR_SPLIT[2])
  print(paste0("t=",t))
}

index_to_choose=which.min(abs(MR_SPLIT_Results[,1] - quantile(MR_SPLIT_Results[,1], 0.5)))
final=MR_SPLIT_Results[index_to_choose,]

#final estimation
final

#index of invalid IVs
G_invalid[[index_to_choose]]

#index of relevant IVs
MR_relevant[[index_to_choose]]

#save results
save(final,MR_SPLIT_Results,MR_relevant,G_invalid,file="MR_SPLIT+.RData")
