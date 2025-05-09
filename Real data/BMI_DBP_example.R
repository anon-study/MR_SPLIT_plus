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

source("../all functions.R")
#all
{
  #gene from GWAS
  plink_file <- "BMI_DBP_10_gene"
  snp_data_G <- snp_attach(paste0(plink_file, ".rds"))
  G_G <- snp_data_G$genotypes
  #39889    21
  
  #gene from samples
  plink_file <- "BMI_DBP_10_all_p001"
  snp_data <- snp_attach(paste0(plink_file, ".rds"))
  G <- snp_data$genotypes
  #39889  4524
  
  #X_Y
  load("BMI_DBP_10_xy.RData")
  
  hist(X)
  hist(log(X))
  
  hist(Y)
  hist(log(Y))
  
  fit=lm(log(Y)~log(X))
  summary(fit)
  #0.224914   0.004037   55.71   <2e-16 
  
  plot(Y~X)
  plot(log(Y)~log(X))
  
  N=nrow(G)
  p=ncol(G)
  
  
  #missing
  missing_ratio <- function(mat) {
    # Check if the input is a matrix
    if (!is.matrix(mat)) stop("Input must be a matrix.")
    
    # calculate the proportion of missing data for each column
    col_missing_ratio <- colSums(is.na(mat)) / nrow(mat)
    
    # return the missing proportions
    return(col_missing_ratio)
  }
  
  miss_G=missing_ratio(G[])
  summary(miss_G)
  #    Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
  #0.0002507 0.0016546 0.0024192 0.0045542 0.0040111 0.0509163 
  
  fill_na_with_random <- function(mat) {
    # Check if the input is a matrix
    if (!is.matrix(mat)) stop("Input must be a matrix.")
    
    # Iterate through each column
    for (col in seq_len(ncol(mat))) {
      # extract the current column
      current_col <- mat[, col]
      
      # calculate the proportion of non-NA values
      value_counts <- table(current_col[!is.na(current_col)])
      proportions <- value_counts / sum(value_counts)
      
      # check if there are 0, 1, or 2; otherwise, set the proportion to 0.
      prop_0 <- ifelse("0" %in% names(proportions), proportions["0"], 0)
      prop_1 <- ifelse("1" %in% names(proportions), proportions["1"], 0)
      prop_2 <- ifelse("2" %in% names(proportions), proportions["2"], 0)
      
      # Generate random values to fill in the NAs
      na_indices <- which(is.na(current_col))
      if (length(na_indices) > 0) {
        mat[na_indices, col] <- sample(
          c(0, 1, 2), 
          size = length(na_indices), 
          replace = TRUE, 
          prob = c(prop_0, prop_1, prop_2)
        )
      }
    }
    
    # return the matrix with imputed values
    return(mat)
  }
  
  set.seed(207)
  G_n=fill_na_with_random(G[])
  
  G_GWAS_n=fill_na_with_random(G_G[])
  
  #save(G_n,G_GWAS_n,X,Y,file="BMI_DBP_10_alldata.RData")
  
  load("BMI_DBP_10_alldata.RData")
  
  #log-transformation
  X=log(X)
  Y=log(Y)
  
  #Center
  X_c=X-mean(X)
  Y_c=Y-mean(Y)
  G_c=apply(G_n, 2, function(col) col - mean(col))
  #columns of G scale to unit length
  G_unit <- apply(G_c, 2, function(col) col / sqrt(sum(col^2)))
  rm(G_c,G_n)
  
  repeat_time=30#repeat times of MR-SPLIT+
  k=2#split the samples into to subsets when applying MR-SPLIT+
  set.seed(927)
  
  #SIS  
  fit=screening(G_unit,X_c,method = "sis",num.select = 80)
  SIS_select=fit$screen
  select_G_sample=snp_data$map[,2][SIS_select]
  select_G_GWAS=snp_data_G$map[,2]
  both_select=which(select_G_sample %in% select_G_GWAS)
  #2 IVs in common
  
  G_candidate=cbind(G_unit[,SIS_select][,-both_select],G_GWAS_n[])
  colnames(G_candidate)=c(select_G_sample[-both_select],select_G_GWAS)
  
  #Standardize gene matrix
  G_c=apply(G_candidate, 2, function(col) col - mean(col))
  #columns of G scale to unit length
  G_candidate <- apply(G_c, 2, function(col) col / sqrt(sum(col^2)))
  rm(G_c,G_unit,G_GWAS_n)
  N_n=nrow(G_candidate)
  p_n=ncol(G_candidate)
  #39889    99
  
  #MR-SPLIT+
  hatX=matrix(0,ncol=repeat_time,nrow=N_n)
  MR_relevant=list()
  selectIV_results=matrix(0,ncol = 6,nrow = repeat_time)
  MR_SPLIT_Results=matrix(0,ncol = 5,nrow = repeat_time)
  colnames(MR_SPLIT_Results)=c("Estimation","Std","pvalue","lc","uc")
  G_invalid=list()#invalid IVs
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
  final
  G_invalid[[index_to_choose]]
  MR_relevant[[index_to_choose]]
  save(final,MR_SPLIT_Results,MR_relevant,G_invalid,file="MR_SPLIT_BMI_DBP_10.RData")
 
  #CIIV
  CI_results = CIIV(Y_c,X_c,G_candidate,robust=FALSE, firststage=TRUE)
  CI_estimation=c(CI_results$Coefficients_CIM[1],CI_results$sd_CIM[1],
                  1-2*abs(pnorm(CI_results$Coefficients_CIM[1],sd=CI_results$sd_CIM[1])-0.5),
                  CI_results$ci_CIM)
  
  ## TSLS
  fit=lm(X_c~G_candidate)
  F_stat=summary(fit)$fstat[1]
  hatX_tsls=fit$fitted.values
  fit2=lm(Y_c~hatX_tsls)
  TSLS=summary(fit2)$coef[2,c(1,2,4)]
  TSLS=c(TSLS,"lc"=TSLS[1]-qnorm(0.975)*TSLS[2],"uc"=TSLS[1]+qnorm(0.975)*TSLS[2])
  names(CI_estimation)=names(TSLS)
  
  #sisVIVE
  sisVIVE_resluts = cv.sisVIVE(Y_c,X_c,G_candidate)
  sisVIVE_estimate = sisVIVE_resluts$beta
  
  #WIT 
  WIT_Results = WIT_practice(X_c,Y_c,G_candidate,seq(0.01,0.2,0.02),ini_lam = 0.05,num_trail = 3)
  WIT_Results_all=WIT_Results[1]$Final
  WIT_Results=c("Est"=WIT_Results_all$WIT,"Std"=WIT_Results_all$WIT_robust_std,
                "pvalue"=1-2*abs(pnorm(WIT_Results_all$WIT,0,WIT_Results_all$WIT_robust_std)-0.5),
                "Lc"=WIT_Results_all$WIT_CI[1],"Uc"=WIT_Results_all$WIT_CI[2])
  
  all_results=list("MR-SPLIT+"=final,
                   "TSLS"=TSLS,
                   "sisVIVE"=sisVIVE_estimate,
                   "CIIV"=CI_estimation,
                   "WIT"=WIT_Results
  )
  all_results
  save(all_results,file="BMI_DBP_10_results.RData")
}


#young
{
  #gene from GWAS
  plink_file <- "BMI_DBP_young_01_GWAS"
  snp_data_G <- snp_attach(paste0(plink_file, ".rds"))
  G_G <- snp_data_G$genotypes
  #12022    21
  
  #gene from samples
  plink_file <- "BMI_DBP_young_01_sample"
  snp_data <- snp_attach(paste0(plink_file, ".rds"))
  G <- snp_data$genotypes
  #12022  3128
  
  #X_Y
  load("BMI_DBP_young_10_xy.RData")
  
  hist(X)
  hist(log(X))
  
  hist(Y)
  hist(log(Y))
  
  fit=lm(log(Y)~log(X))
  summary(fit)
  #0.281395   0.006664   42.23   <2e-16
  
  plot(Y~X)
  plot(log(Y)~log(X))
  
  N=nrow(G)
  p=ncol(G)
  
  
  #missing
  miss_G=missing_ratio(G[])
  summary(miss_G)
  #    Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
  #8.318e-05 1.580e-03 2.412e-03 4.398e-03 3.993e-03 5.116e-02 
  
  set.seed(207)
  G_n=fill_na_with_random(G[])
  G_GWAS_n=fill_na_with_random(G_G[])
  
  #save(G_n,G_GWAS_n,X,Y,file="BMI_DBP_10_young_alldata.RData")
  load("BMI_DBP_10_young_alldata.RData")
  
  #Log-transformation
  X=log(X)
  Y=log(Y)
  #Center
  X_c=X-mean(X)
  Y_c=Y-mean(Y)
  G_c=apply(G_n, 2, function(col) col - mean(col))
  #columns of G scale to unit length
  G_unit <- apply(G_c, 2, function(col) col / sqrt(sum(col^2)))
  rm(G_c,G_n)
  
  set.seed(927)
  #SIS  
  fit=screening(G_unit,X_c,method = "sis",num.select = 80)
  SIS_select=fit$screen
  select_G_sample=snp_data$map[,2][SIS_select]
  select_G_GWAS=snp_data_G$map[,2]
  both_select=which(select_G_sample %in% select_G_GWAS)
  #0 in common
  
  G_candidate=cbind(G_unit[,SIS_select],G_GWAS_n[])
  colnames(G_candidate)=c(select_G_sample,select_G_GWAS)
  G_c=apply(G_candidate, 2, function(col) col - mean(col))
  #columns of G scale to unit length
  G_candidate <- apply(G_c, 2, function(col) col / sqrt(sum(col^2)))
  rm(G_c,G_unit,G_GWAS_n)
  N_n=nrow(G_candidate)
  p_n=ncol(G_candidate)
  #12022 101
  
  #MR-SPLIT+
  repeat_time=30#repeat times of MR-SPLIT+
  k=2#split the samples into to subsets when applying MR-SPLIT+
  #stage one
  hatX=matrix(0,ncol=repeat_time,nrow=N_n)
  MR_relevant=list()
  selectIV_results=matrix(0,ncol = 6,nrow = repeat_time)
  MR_SPLIT_Results=matrix(0,ncol = 5,nrow = repeat_time)
  colnames(MR_SPLIT_Results)=c("Estimation","Std","pvalue","lc","uc")
  G_invalid=list()
  for(t in 1:repeat_time){
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
  final
  G_invalid[[index_to_choose]]
  MR_relevant[[index_to_choose]]
  save(final,MR_SPLIT_Results,MR_relevant,G_invalid,file="MR_SPLIT_BMI_DBP_young_10.RData")
 
  #CIIV
  CI_results = CIIV(Y_c,X_c,G_candidate,robust=FALSE, firststage=TRUE)
  CI_estimation=c(CI_results$Coefficients_CIM[1],CI_results$sd_CIM[1],
                  1-2*abs(pnorm(CI_results$Coefficients_CIM[1],sd=CI_results$sd_CIM[1])-0.5),
                  CI_results$ci_CIM)
 
  ## TSLS
  fit=lm(X_c~G_candidate)
  F_stat=summary(fit)$fstat[1]
  hatX_tsls=fit$fitted.values
  fit2=lm(Y_c~hatX_tsls)
  TSLS=summary(fit2)$coef[2,c(1,2,4)]
  TSLS=c(TSLS,"lc"=TSLS[1]-qnorm(0.975)*TSLS[2],"uc"=TSLS[1]+qnorm(0.975)*TSLS[2])
  names(CI_estimation)=names(TSLS)
  
  #sisVIVE
  sisVIVE_resluts = cv.sisVIVE(Y_c,X_c,G_candidate)
  sisVIVE_estimate = sisVIVE_resluts$beta
 
  #WIT 
  WIT_Results = WIT_practice(X_c,Y_c,G_candidate,seq(0.01,0.2,0.02),ini_lam = 0.05,num_trail = 3)
  WIT_Results_all=WIT_Results[1]$Final
  WIT_Results=c("Est"=WIT_Results_all$WIT,"Std"=WIT_Results_all$WIT_robust_std,
                "pvalue"=1-2*abs(pnorm(WIT_Results_all$WIT,0,WIT_Results_all$WIT_robust_std)-0.5),
                "Lc"=WIT_Results_all$WIT_CI[1],"Uc"=WIT_Results_all$WIT_CI[2])
  
  all_results=list("MR-SPLIT+"=final,
                   "TSLS"=TSLS,
                   "sisVIVE"=sisVIVE_estimate,
                   "CIIV"=CI_estimation,
                   "WIT"=WIT_Results
  )
  all_results
  save(all_results,file="BMI_DBP_10_young_results.RData")
}