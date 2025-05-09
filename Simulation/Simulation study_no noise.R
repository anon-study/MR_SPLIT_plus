rm(list=ls())
#load R packages
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
source("../all functions.R")

#simulation settings
J=21#number of SNPs
N=c(1000, 3000)#sample size
beta=0 #causal effect
J1=21#number of IVs
Times=1000#repeat times of the simulation
repeat_time=10#repeat times of MR-SPLIT+
k=2#split sample into k subsets

#IVs' direct effects on exposure
alpha_1=matrix(c(rep(0.4,21), 
                 rep(0.15,21),
                 c(rep(0.15,5),rep(0.07,16)),
                 rep(0.07,21)),J1,4)

#IVs' direct effects on outcome
alpha_2=matrix(c(c(rep(0,9),rep(0.4,6),rep(0.2,6)),
                 c(rep(0,9),rep(0.4,6),rep(0.2,6)),
                 c(0,0,0,0.2,0.1,rep(0,6),rep(0.2,5),rep(0.1,5)),
                 c(rep(0,9),rep(0.2,6),rep(0.1,6))),J1,4)

#simulations
cl <- makeCluster(10)
registerDoParallel(cl)
#sample size
for(i1 in 1:2){
  #effect size of SNPs
  for(i2 in 1:4){
    results <- foreach(i3=1:Times, .combine=rbind, .packages=c('ivreg','MASS','glmnet','sisVIVE','Matrix','screening','bestsubset')) %dopar% {
      source("../all functions.R")
      
      #Generate samples
      set.seed(i3)
      gamma=alpha_1[,i2]
      Gamma=alpha_2[,i2]
      Sigma=matrix(c(0.5^2,0.3,0.3,1),nrow=2,ncol=2)
      data_set=generate_sample(beta,gamma,Gamma,N[i1],J,Sigma)
      X=data_set$X
      Y=data_set$Y
      G=data_set$G
      real=data_set$real_IV#index of IVs
      
      valid_codes=real[Gamma==0]#index of valid IVs
      invalid_codes=real[Gamma!=0]#index of invalid IVs
      
      #Center the variables
      X_c=X-mean(X)
      Y_c=Y-mean(Y)
      
      #standardize gene matrix
      G_c=apply(G, 2, function(col) col - mean(col))
      #columns of G scale to unit length
      G_unit <- apply(G_c, 2, function(col) col / sqrt(sum(col^2)))
      
      hatX=matrix(0,ncol=repeat_time,nrow=N[i1])
      MR_relevant=list()
      selectIV_results=matrix(0,ncol = 6,nrow = repeat_time)
      MR_SPLIT_Results=matrix(0,ncol = 5,nrow = repeat_time)
      colnames(MR_SPLIT_Results)=c("Estimation","Std","pvalue","lc","uc")
      for(t in 1:repeat_time){
        set.seed(1000+(i3-1)*50+t)
        #stage one of MR-SPLIT+
        MR_stageone=MR_SPLIT_one_nonoise(X_c,Y_c,G_unit,k=2,F_threshold=1000000)
        hatX[,t]=MR_stageone$hatX
        
        MR_relevant[[t]]=c(1:21)
        
        #stage two of MR-SPLIT+
        Z_trans <- (diag(length(hatX[,t])) - hatX[,t] %*% solve(t(hatX[,t]) %*% hatX[,t]) %*% t(hatX[,t])) %*% G_unit[,MR_relevant[[t]]]
        
        #test if there are any invalid IVs  
        many_est = IVreg(Y_c~-1+X_c|G_unit[,MR_relevant[[t]]],  inference = "re") ## Depends on ManyIV package
        test = IVoverid(many_est)
        if(test[2,2]<0.05){
          for(ii in 1:(ncol(Z_trans)-1)){
            results=bs(Z_trans,Y_c,k=ii)
            invalid_index=which(results$beta!=0)
            valid_index=which(results$beta==0)
            many_est = IVreg(Y_c~-1+X_c+G_unit[,MR_relevant[[t]][invalid_index]]|G_unit[,MR_relevant[[t]][valid_index]],  inference = "re") ## Depends on ManyIV package
            test = IVoverid(many_est)
            if(test[2,2]>0.05)break
          }
          optimal_k=ii
          results=bs(Z_trans,Y_c,k=optimal_k)
          
          G_invalid=MR_relevant[[t]][results$beta!=0]
          
        }else{
          G_invalid=NULL
        }
        MR_SPLIT_valid=identify_components_percentage(setdiff(MR_relevant[[t]],G_invalid),valid_codes,invalid_codes)
        MR_SPLIT_invalid=identify_components_percentage(G_invalid,valid_codes,invalid_codes)
        
        selectIV_results[t,]=c(MR_SPLIT_valid, MR_SPLIT_invalid)
        if(length(G_invalid)>0){
          fit=lm(Y_c~hatX[,t]+ G_unit[,G_invalid])
        }else{
          fit=lm(Y_c~hatX[,t])
        }
        MR_SPLIT=c(summary(fit)$coef[2,c(1,2,4)])
        MR_SPLIT_Results[t,]=c(MR_SPLIT,MR_SPLIT[1]-qnorm(0.975)*MR_SPLIT[2],MR_SPLIT[1]+qnorm(0.975)*MR_SPLIT[2])
      }
      results_all=list("MR_SPLIT"=MR_SPLIT_Results,"IV_results"=selectIV_results)
      save(results_all,file=paste0("Results/multiple10_MR_SPLIT+_no_beta0 N_",i1,"_SNP",i2,"_Times",i3,".RData"))
    }
  }
}  
stopCluster(cl)