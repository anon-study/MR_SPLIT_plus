## WIT-Estimator 
library(ncvreg)
## Use the following to install
#install.packages("ncvreg")
library(ManyIV)
## Use the following to install 
# devtools::install_github("https://github.com/kolesarm/ManyIV")
library(ivmodel)
# use the following to install 
# install.packages("ivmodel")
#########################################################
##################### Description #######################
#########################################################
# In the following description
# Y is a nx1 vector standards for responds
# Z is a nxp vector means potential p-IVs that contains 
# Valid and invalid IVs
# D is a nx1 vector means endogenous variables of interest
#########################################################
## Helper function
## MCP with two parameter lambda and kappa, kappa = 1/alpha, which alpha is usually = 2;
mcp = function(x,lambda,kappa)
{
  (lambda*abs(x)-kappa*x^2/2)*(abs(x)<=lambda/kappa) + (lambda^2/(2*kappa))*(abs(x)>lambda/kappa)
}
## Derivative of MCP penalty 
dmcp = function(x,lambda,kappa)
{
  (lambda-kappa*abs(x))*(abs(x)<=(lambda/kappa))
}
## Operator of SoftThresholding
SoftThresholding = function(X,Lambda)
{
  sign(X)*pmax(0,abs(X)-Lambda)
}
dL = function(nXX,nXY,beta0)
{
  (-nXY+nXX%*%beta0)
}
## Transformed Loss
L = function(nXX,nXY,nYY,x)
{
  nYY-2*crossprod(x,nXY)+crossprod(x,(nXX%*%x))
}
FPR_FNR = function(true_valid_index,estimate_valid)
{
  p = length(estimate_valid)
  true_invalid_index = rep(1,p)-true_valid_index
  estimate_invalid = rep(1,p)-estimate_valid
  
  TP = sum(true_valid_index*estimate_valid)
  FP = sum(true_invalid_index*estimate_valid)
  FN = sum(true_valid_index*estimate_invalid)
  TN = sum(true_invalid_index*estimate_invalid)
  
  result = list("FPR" = (FP/(FP+TN)),"FNR" = (FN/(FN+TP)))
  return(result)
}
ALL_oracle = function(alpha,PI)
{
  uni = unique(alpha/PI)
  number_oracle = length(uni)
  All_oracle_solutions = matrix(rep(0,length(alpha)*number_oracle),number_oracle,length(alpha))
  for (i in 1:number_oracle)
  {
    All_oracle_solutions[i,] = t(alpha+PI*-uni[i])
  }
  return(All_oracle_solutions)
}
######################
## I-LAMM algorithm ##
######################
## LAMM updating (inner iteration)
LAMM = function(x_old,CAP_Lambda,phi,kappa,nYY,nXX,nXY)
{
  x_new = SoftThresholding(x_old-dL(nXX,nXY,x_old)/phi,CAP_Lambda/phi)
  return(x_new)
}
## I-LAMM updating (outer iteration)
ILAMM = function(CAP_Lambda,x0,phi,kappa,nYY,nXX,nXY,error)
{
  x_old = x0
  TMAX = 500;
  for (J in 1:TMAX)
  {
    x_new = LAMM(x_old,CAP_Lambda,phi,kappa,nYY,nXX,nXY)
    judge = (SoftThresholding(dL(nXX,nXY,x_new),abs(CAP_Lambda)))*(x_new ==0)+(dL(nXX,nXY,x_new)+CAP_Lambda*sign(x_new))*(x_new!=0)
    
    if(norm(judge,"I")<= error)
    {
      print(paste("sub-iteration = ", J))
      
      return(x_new)
      break
    }
    x_old = x_new
  }
  return(x_new)
}
## FUll-ILAMM algorithm that integrates the previous two terms
FINAL_ILAMM = function(Y,X,beta_old,lambda,kappa,error_c,error_t)
{
  
  p = ncol(X)
  n = nrow(X)
  nYY = t(Y)%*%Y/n
  nXY = t(X)%*%Y/n
  nXX = t(X)%*%X/n
  iterationMAX = 2000
  phi = max(eigen(nXX)$value)## only compute once to achieve the majorization.
  ## Outer iteration
  for(i in 1: iterationMAX)
  {
    print(paste("iteration = ", i))
    CAP_Lambda = dmcp(abs(beta_old),lambda,kappa)
    if (i==1)
    {
      error = error_c
    }
    else
    {
      error = error_t
    }
    ## complete
    beta_new = ILAMM(CAP_Lambda,beta_old,phi,kappa,nYY,nXX,nXY,error)
    
    if(norm(beta_new-beta_old,"I") < error )
    {
      print("success I-LAMM algorithm")
      Loss = L(nXX,nXY,nYY,beta_new)
      print(paste("Loss = ", Loss))
      return(beta_new)
      break
    }
    beta_old = beta_new
  }
  return(beta_new)
}
######################
## Main Body of WIT ##
######################
## WIT estimator with a predetermined initial point and a fixed tuning parameter lambda
WIT = function(D,Y,Z,beta_old,lambda,kappa=0.5,error_c=1e-3,error_t=1e-5)
{
  if( (!is.vector(Y) && !is.matrix(Y)) | !is.numeric(Y)) stop("Y is not a numeric vector.")
  if( (!is.vector(D) && !is.matrix(D)) | !is.numeric(D)) stop("D is not a numeric vector.")
  if( !is.matrix(Z) | !is.numeric(Z)) stop("Z is not a numerical matrix.")
  if(nrow(Z) != length(Y)) stop("The dimension of Y differs from the row dimension of Z")
  if(nrow(Z) != length(D)) stop("The dimension of D differs from the row dimension of Z")
  if(length(D) != length(Y)) stop("The dimension of Y differs from the dimension of D")
  
  n = nrow(Z); L = ncol(Z);
  p = L
  
  ## Compute transformed tilde_Z and tilde_Y
  Yhat = qr.fitted(qr(Z),Y); Dhat = qr.fitted(qr(Z),D);
  
  Y.transformed = Yhat - Dhat * (sum(Dhat * Yhat) / sum( Dhat^2))
  
  Z.transformed = qr.resid(qr(Dhat),Z)#Z - Dhat %*% t(Dhat) %*% Z / sum( Dhat^2)
  
  ### SELECTION0-STAGE
  ## Using I-LAMM algorithm to obtain alpha_MCP  
  alpha_hat = FINAL_ILAMM(Y.transformed,Z.transformed,beta_old,lambda,kappa,error_c,error_t)
  rm(Z.transformed)
  rm(Y.transformed)
  
  if (sum(alpha_hat !=0)==p)
  {
    return(NA)
  }
  alphaSuppSize = apply(alpha_hat,1,function(x){sum(x != 0)}); indexProperSupp = which(alphaSuppSize < p); 
  
  invalid_index = which(alphaSuppSize ==1)
  valid_index = which(alphaSuppSize ==0)
  ##
  
  Z_valid = Z[,valid_index]
  Z_invalid = Z[,invalid_index]
  
  # judge for if no zeros: 
  if (length(valid_index) ==0)
  {
    print("Unidentified")
  }else
  {
    
    ### ESTIMATION-STAGE using LIML
    if (length(invalid_index)!=0)
    {
      
      many_est = IVreg(Y~-1+D+Z[,invalid_index]|Z[,valid_index],  inference = "re") ## Depends on ManyIV package
      xxx = ivmodel(Y,D,Z[,valid_index],Z[,invalid_index],intercept=FALSE) # for LIML alpha
      WIT =  many_est$estimate[3,1] #as.numeric(xxx$Fuller$point.est)      #many_est$estimate[3,1] ## Liml 
      WIT_robost_std = many_est$estimate[3,2]
      test = IVoverid(many_est)
      
      ### Get LIML result of alpha estimate
      WIT_alpha = rep(0,p); 
      WIT_alpha[invalid_index] = xxx$LIML$point.est.other
      
      WIT_CI = c(WIT-1.96*WIT_robost_std,WIT+1.96*WIT_robost_std)
      ## Summary of Result
      Result = list("WIT" = WIT,"WIT_robust_std" = WIT_robost_std,"WIT_CI" = WIT_CI,"alpha_MCP" = as.numeric(alpha_hat),"WIT_valid" = valid_index,"alpha_WIT" = WIT_alpha,"Sargen_p_value" = test[1,2],"Modified-CD"= test[2,2])
      return(Result)
    }
    else
    {
      many_est = IVreg(Y~-1+D|Z[,valid_index],  inference = "re") ## Depends on ManyIV package
      xxx = ivmodel(Y,D,Z[,valid_index],intercept=FALSE) # for LIML alpha
      WIT =  many_est$estimate[3,1] #as.numeric(xxx$Fuller$point.est)      #many_est$estimate[3,1] ## liml 
      WIT_robost_std = many_est$estimate[3,2]
      test = IVoverid(many_est)
      
      ### Get LIML result of alpha estimate
      WIT_alpha = rep(0,p); 
      WIT_alpha[invalid_index] = NaN
      
      WIT_CI = c(WIT-1.96*WIT_robost_std,WIT+1.96*WIT_robost_std)
      ## Summary of Result
      Result = list("WIT" = WIT,"WIT_robust_std" = WIT_robost_std,"WIT_CI" = WIT_CI,"alpha_MCP" = as.numeric(alpha_hat),"WIT_valid" = valid_index,"alpha_WIT" = WIT_alpha,"Sargen_p_value" = test[1,2],"Modified-CD"= test[2,2])
      return(Result)
    }
  }
}
######################################
## Main Body of WIT with MCD tuning ##
######################################
## Tuning Strategy with MCD Test
MCD_tuning = function(D,Y,Z,Lambda_Seq,kappa,error_c,error_t,initial) ## 
{ 
  
  n=dim(Z)[1]
  ## Record MCD value
  MCD_value = rep(0,length(Lambda_Seq))
  Model_size_value = rep(0,length(Lambda_Seq))
  
  ## For Record
  History = list()
  for (ii in (1:length(Lambda_Seq)))
  {
    lambda = Lambda_Seq[ii]
    
    WIT_Current = WIT(D,Y,Z,initial,lambda,kappa,error_c,error_t)
    History[[length(History)+1]] = WIT_Current
    if (all(is.na(WIT_Current)))
    {
      next
    }
    Model_size_value[ii] = sum(WIT_Current$alpha_WIT !=0)
    MCD_value[ii] = WIT_Current$`Modified-CD`
  }
  ## Selected the result 
  filtered_index = (MCD_value > 0.5/log(n))
  if (sum(is.na(MCD_value)) == length(Lambda_Seq)){
    # not censored
    filtered_index = 1:length(Lambda_Seq)
  }
  
  filtered_index[is.na(filtered_index)] = rep(T,sum(is.na(filtered_index)))
  
  if(sum(!filtered_index) == length(Lambda_Seq)) ## all cannot pass MCD
  {
    return(NA)
  }
  
  Optimal_Model_size_value = which(Model_size_value == min(Model_size_value[filtered_index]))[1]
  
  WIT_final_estimation = WIT(D,Y,Z,initial,Lambda_Seq[Optimal_Model_size_value],kappa,error_c,error_t)
  return(list("Chosen_Model" = WIT_final_estimation,"HISTORY" = History))
}
WIT_practice = function(D,Y,Z,Lambda_Seq,ini_lam = 0.05,kappa=2,error_c=1e-3,error_t=1e-4,num_trail=2)
{
  ## Ratio Estimate
  n= dim(Z)[1]
  p = dim(Z)[2]
  ## can be modified by ini_lam
  Gamma_est = qr.solve(Z,Y)
  gamma_est = qr.solve(Z,D)
  
  ## Record History
  
  inidividual_esti = Gamma_est/gamma_est
  ord_ini_est = order (inidividual_esti)
  sorted_ini_est = sort(inidividual_esti)
  XXX = matrix(rep(1,p*p),p,p);
  XXX[lower.tri(XXX)] = 0;
  
  ## Rough clustering
  fuzed_model = ncvfit(XXX,sorted_ini_est,sorted_ini_est,penalty = "MCP",lambda = ini_lam,penalty.factor = c(0,rep(1,length(inidividual_esti)-1)))
  
  fuzed_ini_est = XXX%*%fuzed_model$beta
  rm(fuzed_model)
  uni_fuzed_ini_est = unique(fuzed_ini_est)
  number_group = length(uni_fuzed_ini_est)
  group_count = matrix(c(rep(0,2*number_group)),number_group,2)
  group_count[,1] = uni_fuzed_ini_est
  for (i in 1:length(uni_fuzed_ini_est))
  {
    group_count[i,2] = sum(fuzed_ini_est == uni_fuzed_ini_est[i])
    
  }
  
  ## Determine how many initials we use
  num_trail = min(num_trail,number_group)
  order_in = order(group_count[,2],decreasing = T)
  group_count = group_count[order_in,]
  
  ## Compute results
  initials = matrix(rep(0,p*(num_trail+1)),p,num_trail+1)
  for (j in 1:num_trail)
  {
    initials[,j] = Gamma_est-group_count[j,1]*gamma_est
    initials[(ord_ini_est[fuzed_ini_est==group_count[j,1]]),j] = rep(0,group_count[j,2])
  }
  l_0_alpha = p+1
  
  ## Tuning with MCD
  
  ## Record History
  history = list()
  final_model=NA
  for (j in 1:(num_trail+1))
  { 
    temp = MCD_tuning(D,Y,Z,Lambda_Seq,kappa,error_c,error_t,initials[,j])
    if(length(temp) == 1 && is.na(temp)){
      next
    }
    med_model = temp[1]$Chosen_Model
    history[[length(history)+1]] = temp[2]$HISTORY
    if (all(is.na(med_model)))
    {
      next
    }
    if (sum(med_model$alpha_WIT!=0)<l_0_alpha)
    {
      l_0_alpha = sum(med_model$alpha_WIT!=0)
      final_model = med_model
    }
  }
  return(list("Final" = final_model,"Tuning_History" = history))
}


#CIIV
CIIV <- function(Y, D, Z, X, alpha = 0.05, tuning = 0.1/log(length(Y)), robust = TRUE, firststage = FALSE, firsttuning = sqrt(2.01*log(ncol(Z)))) {
  #Check Input
  #Check Data
  
  if(is.data.frame(Y)){Y <- as.matrix(Y)}
  if(!is.vector(Y) && !is.matrix(Y) | !is.numeric(Y) | ncol(as.matrix(Y))!=1)
    stop("Y must be a numeric vector.");
  Y = as.numeric(Y);
  
  if(!is.null(colnames(D))){Dname = colnames(D)}else{Dname = "D"};
  if(is.data.frame(D)){D <- as.matrix(D)}
  if(!is.vector(D) && !is.matrix(D) | !is.numeric(D) | ncol(as.matrix(D))!=1)
    stop("D must be a numeric vector.");
  D = as.numeric(D);
  
  if(is.data.frame(Z)){Z <- as.matrix(Z)}
  if(!is.matrix(Z) | !is.numeric(Z))
    stop("Z must be a numerical matrix.");
  if(ncol(Z)<2)
    stop("The number of instruments must be greater than 1.");
  if( ncol(Z) > nrow(Z))
    stop("The number of instruments must be smaller than the sample size.");
  
  if(!missing(X) && is.data.frame(X)){X <- as.matrix(X)}
  if(!missing(X) && (!is.matrix(X) | !is.numeric(X)))
    stop("X must be a numerical matrix.");
  
  stopifnot(length(Y) == length(D),length(Y) == nrow(Z));
  
  #Other Arguments
  stopifnot(is.logical(firststage), is.logical(robust));
  stopifnot(is.numeric(alpha), length(alpha) == 1,alpha <= 1,alpha >= 0);
  stopifnot(is.numeric(tuning), length(tuning) == 1);
  
  # Define Constants
  n <- length(Y);
  pz <- ncol(Z);
  
  # Preserve Data
  d = D;
  y = Y;
  z = Z;
  
  if (!missing(X)){
    X <- cbind(1, X);
    if(is.null(colnames(X))){colnames(X) = c("intercept", paste0('X', 1:ncol(X)))};
  }else{
    X <- matrix(1, n, 1);
    colnames(X) <- "intercept";
  };
  Covariates = cbind(D,X);
  Exogenous = cbind(Z,X);
  
  # Variable Names
  colnames(Covariates)[1] <- Dname;
  CovariatesNames <- colnames(Covariates);
  if(!is.null(colnames(Z))){InstrumentNames = colnames(Z)}else{InstrumentNames = paste0('Z', 1:ncol(Z))};
  
  # Centralization
  D <- qr.resid(qr(X), D);
  Z <- scale(qr.resid(qr(X), Z), center = FALSE, scale = TRUE);
  
  # First Stage
  if (firststage){
    
    lm_first <- lm(D ~ Z - 1);
    coef_first <- coef(lm_first);
    sd_first <- sqrt(diag(if(robust) vcovHC(lm_first, "HC0") else vcov(lm_first)));
    
    first_index <- as.numeric(abs(coef_first/sd_first) >= firsttuning);
    relevant_index <- relevant_instruments <- which(first_index == 1);
    Nr_relevant <- length(relevant_instruments);
    
    if (sum(first_index) <= 1){
      warning("Less than two IVs are individually relevant, treat all IVs as strong");
      relevant_instruments <- c(1:pz);
    }
    
    if (length(relevant_instruments) < pz){
      X <- cbind(X, z[, -relevant_instruments]);
      Z <- as.matrix(z[, relevant_instruments], nrow = n, ncol = Nr_relevant);
      pz <- ncol(Z);
      
      D <- qr.resid(qr(X), d);
      Z <- scale(qr.resid(qr(X), Z), center = FALSE, scale = TRUE);
    }
    
  }else{
    relevant_index <- relevant_instruments <- c(1:pz);
    Nr_relevant <- pz;
  };
  
  Y <- qr.resid(qr(X), Y);
  
  # OLS Estimation for the IV-Specific Estimates and Standard Error
  lm_Reduced <- lm(cbind(Y, D) ~ Z - 1);
  
  # IV-Specific Estimates
  gamma_Y <- coef(lm_Reduced)[, 1];
  gamma_D <- coef(lm_Reduced)[, 2];
  betaIV <- gamma_Y / gamma_D;
  
  # Variances by Delta Method
  U <- solve(crossprod(Z));
  Jacobian_betaIV <- cbind(
    diag(c(1 / gamma_D)), diag(c(-gamma_Y / gamma_D^2))
  );
  Covmat_gamma <- if(robust) vcovHC(lm_Reduced, "HC0") else vcov(lm_Reduced);
  sdIV <- sqrt(diag(Jacobian_betaIV %*% Covmat_gamma %*% t(Jacobian_betaIV)));
  
  # CI Downward Testing Sargen/Hansen-J Procedure and Post-Selection Estimation
  CIIV.TestingSelectionEstimation(
    Y, D, Z, U,
    Covariates,
    Exogenous,
    y, z,
    CovariatesNames,
    InstrumentNames,
    betaIV = betaIV,
    sdIV = sdIV,
    alpha = alpha,
    tuning = tuning,
    robust = robust,
    gamma_D = gamma_D,
    Covmat_gamma = Covmat_gamma,
    firststage = firststage,
    relevant_instruments = relevant_instruments,
    Nr_relevant = Nr_relevant,
    relevant_index = relevant_index
  );
}
CIIV.TestingSelectionEstimation <- function(
    Y, D, Z, U, Covariates, Exogenous, y, z, CovariatesNames, InstrumentNames, betaIV, sdIV, alpha = 0.05, tuning = 0.1/log(length(Y)), robust = TRUE,
    gamma_D, Covmat_gamma, firststage = FALSE, relevant_instruments, Nr_relevant, relevant_index) {
  
  # Function for Sargan Test
  non_robust_sar_CI <- function(res_CI, Z, U, n) {
    (t(res_CI) %*% Z %*% U %*% t(Z) %*% res_CI) /
      (t(res_CI) %*% res_CI / n)
  }
  
  # Function for Two Step GMM and Hansen-J Test
  CIM.HansenJTest <- function(Y,X,Z){
    res_FirstStep <- residuals(AER::ivreg(Y ~ X - 1 | Z));
    
    Weight_SecondStep <- crossprod(res_FirstStep * Z);
    
    Coef_SecondStep <- solve(
      t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
    ) %*% t(X) %*% Z %*% solve(Weight_SecondStep) %*% t(Z) %*% Y;
    
    res_SecondStep <- as.vector(Y - X %*% Coef_SecondStep);
    
    sd_SecondStep <- sqrt(diag(solve(
      t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
    ) %*% t(X) %*% Z %*% solve(Weight_SecondStep)%*%crossprod(res_SecondStep * Z)%*%t(
      solve(
        t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
      ) %*% t(X) %*% Z %*% solve(Weight_SecondStep)
    )));
    
    HansenJ_Stat <- t(res_SecondStep) %*% Z %*% solve(Weight_SecondStep) %*%
      t(Z) %*% res_SecondStep;
    
    list(HansenJ_Stat, Coef_SecondStep,sd_SecondStep)
  }
  
  # Define Constants
  n <- length(Y);
  pz <- ncol(Z);
  
  colnames(Z) = paste0(relevant_instruments);
  colnames(z) = paste0(c(1:ncol(z)));
  
  # Break Points
  crit <- matrix(0, pz, pz);
  rownames(crit) = colnames(crit) <- paste0(relevant_instruments);
  for (i in 1:pz) {
    for (j in 1:pz) {
      crit[i, j] <- abs(betaIV[i] - betaIV[j]) / (sdIV[i] + sdIV[j]);
    }
  }
  
  # The Downward Testing Procedure
  maxs <- pz;
  Psar <- 0;
  epsi <- 10^(-7);
  Psi <- max(crit) + epsi;
  
  while (maxs > 0 && Psar < tuning) {
    # The confidence intervals
    ciIV <- matrix(0, pz, 3);
    ciIV[, 1] <- betaIV - Psi * sdIV;
    ciIV[, 2] <- betaIV + Psi * sdIV;
    ciIV[, 3] <- relevant_instruments;
    
    # Ordering the Confidence Intervals by the Lower Ends
    ciIV <- ciIV[order(ciIV[, 1]), ];
    
    # Check Overlap
    grid_CI <- matrix(0, pz, pz)
    for (i in 2:pz) {
      for (j in 1:i) {
        grid_CI[i, j] <- as.numeric(ciIV[i, 1] <= ciIV[j, 2]);
      }
    }
    
    sumCI <- apply(grid_CI, 1, sum);
    selCI <- as.matrix(grid_CI[which(sumCI == max(sumCI)), ]);
    maxs <- max(sumCI);
    
    #Selection
    if(maxs >= 2 && length(selCI) < pz + 1) {  # No tie
      selCI <- cbind(ciIV[, 3], selCI);
      selVec <- sort(selCI[selCI[, 2] == 1, 1]);
      
      Zv <- Z[,!colnames(Z) %in% paste0(selVec)];
      if(robust) {
        sar_CI <- CIM.HansenJTest(Y, cbind(D, Zv), Z)[[1]];
      } else {
        res_CI <- resid(AER::ivreg(Y ~ cbind(D, Zv) - 1 | Z));
        sar_CI <- non_robust_sar_CI(res_CI, Z, U, n);
      }
      Psi <- max(crit[paste0(selVec), paste0(selVec)]) - epsi; # updated psi
      Psar <- pchisq(sar_CI, maxs - 1, lower.tail = FALSE);
      
    } else if(maxs >= 2) {  # Tie
      selCI <- cbind(ciIV[, 3], t(selCI));
      selVec <- matrix(0, maxs, ncol(selCI) - 1);
      Psar <- Psi <- rep(0, ncol(selCI) - 1);
      
      for (i in 1:(ncol(selCI) - 1)) {
        selVec[, i] <- sort(selCI[selCI[, i+1] == 1, 1]);
        
        if(robust){
          sar_CI <- CIM.HansenJTest(Y, cbind(D, Z[,!colnames(Z) %in% paste0(selVec[, i])]), Z)[[1]];
        } else{
          res_CI <- resid(AER::ivreg(Y ~ D + Z[,!colnames(Z) %in% paste0(selVec[, i])] - 1 | Z));
          sar_CI <- non_robust_sar_CI(res_CI, Z, U, n);
        }
        
        Psar[i] <- pchisq(sar_CI, maxs - 1, lower.tail = FALSE);
        Psi[i] <- max(crit[paste0(selVec[, i]), paste0(selVec[, i])]);
      }
      
      index.tie <- match(max(Psar), Psar);
      selVec <- selVec[,index.tie];
      Psar <- Psar[index.tie];
      Psi <- min(Psi) - epsi;  # updated psi
    } else {  # None Valid
      maxs <- 0;
      Psar <- 0;
      selVec <- NULL;
    }
  }
  
  Nr_valid <- length(selVec);
  
  # Post-Selection Estimation
  if (Nr_valid == 0) {
    print("None of the instruments is selected as valid, do OLS.")
    
    regressor_CIM_temp = regressor_CIM <- cbind(Covariates, z);
    length_regressor <- ncol(regressor_CIM)-ncol(z);
    Coefficients_CIM <- qr.coef(qr(regressor_CIM), y)[1:length_regressor];
    res_CIM <- qr.resid(qr(regressor_CIM), y);
  } else {
    # At least one of the IVs is selected as valid.
    
    z_invalid <- matrix(z[,!colnames(z) %in% paste0(selVec)], ncol = (ncol(z) - Nr_valid), nrow = n);
    
    regressor_CIM_temp <- cbind(Covariates, z_invalid);
    regressor_CIM <- cbind(fitted(lm(Covariates[,1] ~ Exogenous)), Covariates[,-1], z_invalid)
    length_regressor <- ncol(regressor_CIM_temp) - ncol(z_invalid);
    iv_CIM <- AER::ivreg(y ~ regressor_CIM_temp - 1 | Exogenous);
    Coefficients_CIM <- coef(iv_CIM)[1:length_regressor];
    if(robust){
      Coefficients_CIM_GMM <- (CIM.HansenJTest(y, regressor_CIM_temp, Exogenous)[[2]])[1:length_regressor];
    }
    res_CIM <- resid(iv_CIM);
  }
  
  # Standard Error and Confidence Interval
  if (robust) {
    sd_CIM <- sqrt(diag(
      solve(crossprod(regressor_CIM)) %*%
        crossprod(res_CIM * regressor_CIM) %*%
        solve(crossprod(regressor_CIM))
    ))[1:length_regressor];
    ci_CIM <- c(
      Coefficients_CIM[1] - qnorm(1-alpha/2) * sd_CIM[1],
      Coefficients_CIM[1] + qnorm(1-alpha/2) * sd_CIM[1]
    );
    if (Nr_valid == 0){
      Coefficients_CIM_GMM <- NA;
      sd_CIM_GMM <- NA;
      ci_CIM_GMM <- NA;
    }else{
      sd_CIM_GMM <- (CIM.HansenJTest(y, regressor_CIM_temp, Exogenous)[[3]])[1:length_regressor];
      ci_CIM_GMM <- c(
        Coefficients_CIM_GMM[1] - qnorm(1-alpha/2) * sd_CIM_GMM[1],
        Coefficients_CIM_GMM[1] + qnorm(1-alpha/2) * sd_CIM_GMM[1]
      );
    }
  } else {
    sd_CIM <- sqrt(diag(
      mean(res_CIM^2) * solve(crossprod(regressor_CIM))
    ))[1:length_regressor];
    ci_CIM <- c(
      Coefficients_CIM[1] - qnorm(1-alpha/2) * sd_CIM[1],
      Coefficients_CIM[1] + qnorm(1-alpha/2) * sd_CIM[1]
    );
  }
  
  # Results
  if(robust){
    object <- list(
      robust = robust,
      if(is.null(selVec)){
        Valid_Instruments = paste0("None instruments selected as valid.");
      }else{
        Valid_Instruments = InstrumentNames[c(1:length(z))[selVec]];
      },
      Valid_Instruments = Valid_Instruments,
      Nr_valid = Nr_valid,
      if(firststage){
        if(Nr_relevant == 0){
          Relevant_Instruments = paste0("None instruments selected as relevant.")
        }else{
          Relevant_Instruments = InstrumentNames[c(1:length(z))[relevant_index]];
        };
      }else{
        Relevant_Instruments = paste0("No first stage selection.");
        Nr_relevant = paste0("No first stage selection.");
      },
      Relevant_Instruments = Relevant_Instruments,
      Nr_relevant = Nr_relevant,
      Covariate_Names = CovariatesNames,
      Coefficients_CIM = Coefficients_CIM,
      sd_CIM = sd_CIM,
      Coefficients_CIM_GMM = Coefficients_CIM_GMM,
      sd_CIM_GMM = sd_CIM_GMM,
      ci_CIM = ci_CIM,
      ci_CIM_GMM = ci_CIM_GMM,
      HansenJ_CIM = Psar
    )
  }else{
    object <- list(
      robust = robust,
      if(is.null(selVec)){
        Valid_Instruments = paste0("None instruments selected as valid.");
      }else{
        Valid_Instruments = InstrumentNames[c(1:ncol(z))[selVec]];
      },
      Valid_Instruments = Valid_Instruments,
      Nr_valid = Nr_valid,
      if(firststage){
        if(Nr_relevant == 0){
          Relevant_Instruments = paste0("None instruments selected as relevant.")
        }else{
          Relevant_Instruments = InstrumentNames[c(1:length(z))[relevant_index]];
        };
      }else{
        Relevant_Instruments = paste0("No first stage selection.");
        Nr_relevant = paste0("No first stage selection.");
      },
      Relevant_Instruments = Relevant_Instruments,
      Nr_relevant = Nr_relevant,
      Covariate_Names = CovariatesNames,
      Coefficients_CIM = Coefficients_CIM,
      sd_CIM = sd_CIM,
      ci_CIM = ci_CIM,
      Sargan_CIM = Psar
    )
  }
  
  class(object) <- "CIIV"
  
  object
}
print.CIIV <- function(object,robust = object$robust,...){
  cat("\nValid Instruments:\n", object$Valid_Instruments, "\n","\nNumber of Valid Instruments:\n", object$Nr_valid,"\n");
  cat("\nRelevant Instruments:\n", object$Relevant_Instruments, "\n","\nNumber of Relevant Instruments:\n", object$Nr_relevant,"\n");
  cat("\nCoefficients:\n")
  
  names(object$Coefficients_CIM) = names(object$sd_CIM) = names(object$ci_CIM) = NULL;
  if(robust){
    names(object$Coefficients_CIM_GMM) = names(object$sd_CIM_GMM) = names(object$ci_CIM_GMM) = NULL;
    coef_cim <- cbind(object$Coefficients_CIM, object$sd_CIM, object$Coefficients_CIM_GMM,object$sd_CIM_GMM);
    colnames(coef_cim) = c("2SLS Estimate", "2SLS Std. Error", "GMM Estimate", "GMM Std. Error");
    rownames(coef_cim) = object$Covariate_Names;
    print(coef_cim, quote = FALSE);
    
    cat("\nConfidence Interval 2SLS: [", object$ci_CIM[1], ",", object$ci_CIM[2], "]", "\n", sep = '');
    cat("\nConfidence Interval GMM: [", object$ci_CIM_GMM[1], ",", object$ci_CIM_GMM[2], "]", "\n", sep = '');
    
    cat("\np-value of Hansen-J: ", object$HansenJ_CIM, sep = '');
  }else{
    coef_cim <- cbind(object$Coefficients_CIM, object$sd_CIM);
    colnames(coef_cim) = c("2SLS Estimate", "Std. Error");
    rownames(coef_cim) = object$Covariate_Names;
    print(coef_cim, quote = FALSE);
    
    cat("\nConfidence Interval: [", object$ci_CIM[1], ",", object$ci_CIM[2], "]", "\n", sep = '');
    
    cat("\np-value of Sargan: ", object$Sargan_CIM, sep = '');
  }
}
cross_validation <- function(G, Y, K = 10, seed = 123) {
  # 确保输入数据为矩阵和向量
  G <- as.matrix(G)
  Y <- as.vector(Y)
  
  # 获取样本数量和最大 k 值
  n <- nrow(G)
  k_max <- ncol(G)
  
  # 设置随机种子以保证结果可重复
  set.seed(seed)
  
  # 将样本随机分成 K 份
  folds <- sample(rep(1:K, length.out = n))
  
  # 初始化存储残差平方和的向量
  cv_errors <- numeric(k_max)
  cv_errors_sd <- numeric(k_max)
  # 交叉验证循环
  for (k in 1:k_max) {
    fold_errors <- numeric(K)
    
    for (i in 1:K) {
      # 划分训练集和验证集
      train_idx <- which(folds != i)
      valid_idx <- which(folds == i)
      
      G_train <- G[train_idx, ]
      Y_train <- Y[train_idx]
      
      G_valid <- G[valid_idx, ]
      Y_valid <- Y[valid_idx]
      
      # 在训练集上估计 beta
      results <- bs(G_train, Y_train,k)
      beta_k <- results$beta
      
      # 在验证集上计算残差平方和
      residuals <- Y_valid - G_valid %*% beta_k
      fold_errors[i] <- sum(residuals^2)/length(residuals)
    }
    
    # 计算 k 值对应的平均残差平方和
    cv_errors[k] <- mean(fold_errors)
    cv_errors_sd[k]<- sd(fold_errors)
  }
  
  # 返回每个 k 值对应的残差平方和
  return(cbind(1:k_max,cv_errors,cv_errors_sd))
}

#MR_SPLIT+
#Include noise/Real data
MR_SPLIT_one_CIfirst=function(X,Y,G,k=2,F_threshold=1000000){
  N=length(X)
  #split sample
  selected=NULL
  unselected=1:N
  J=ncol(G)
  n=sample(unselected,size=N/k,replace = FALSE)
  selected=cbind(selected,n)
  for(p in 1:k){
    if(p<k){
      n=sample(unselected[-selected],size=N/k,replace = FALSE)
      selected=cbind(selected,n)
    }
  }
  IV_matrix=matrix(0,ncol=4,nrow=k)
  
  selectIV=list()
  weights=list()
  
  for(p in 1:k){
    n=selected[,p]
    
    #LASSO
    #select IVs 
    lambda=cv.glmnet(G[-n,],X[-n],alpha=1)$lambda.min
    lasso=glmnet(G[-n,],X[-n],lambda=lambda)
    #coef of all IVs
    a=coef(lasso)[-1,]
    while(sum(a!=0)==0){
      lambda=0.95*lambda
      lasso=glmnet(G[-n,],X[-n],lambda=lambda)
      a=coef(lasso)[-1,]
    }
    #selected IVs
    selectIV[[p]]=which(a!=0)
    
    #CIIV_firststage
    firsttuning = sqrt(2.01*log(ncol(G[-n,selectIV[[p]],drop=FALSE])))
    #select IVs 
    lm_first <- lm(X[-n] ~ G[-n,selectIV[[p]]] - 1)
    coef_first <- coef(lm_first)
    sd_first <- sqrt(diag(vcov(lm_first)))
    
    first_index <- as.numeric(abs(coef_first/sd_first) >= firsttuning)
    if (sum(first_index) <= 1){
      warning("Less than two IVs are individually relevant, treat all IVs as strong")
      selectIV[[p]] <- selectIV[[p]]
    }else{
      selectIV[[p]] <-  selectIV[[p]][which(first_index == 1)]
    }
    
    fit=lm(X[-n]~G[-n,selectIV[[p]]])
    weights[[p]]=coef(fit)[-1]
    IV_matrix[p,3]=length(selectIV[[p]])
  }
  
  
  #select major and weak based on partial F statistics
  partialF=list()
  for(p in 1:k){
    n=selected[,p]
    full=lm(X[-n]~G[-n,selectIV[[p]]])
    Fs=NULL
    if(length(selectIV[[p]])==1){
      Fs=summary(lm(X[-n]~G[-n,selectIV[[p]]]))$fstat[1]
    }
    if(length(selectIV[[p]])>1){
      for(q in 1:length(selectIV[[p]])){
        reduced=lm(X[-n]~G[-n,selectIV[[p]][-q]])
        Fs[q]=anova(reduced,full)$F[2]
      }
    }
    partialF[[p]]=Fs
  }
  
  selectmajor=list()
  selectweak=list()
  for(p in 1:k){
    if(is.na(selectIV[p])==FALSE){
      selectmajor[[p]]=selectIV[[p]][partialF[[p]]>=F_threshold]
      selectweak[[p]]=selectIV[[p]][partialF[[p]]<F_threshold]
      IV_matrix[p,1]=sum(partialF[[p]]>=F_threshold)
      IV_matrix[p,2]=sum(partialF[[p]]<F_threshold)
    }
  }
  
  #Get coef
  combineMajor=list()
  for(p in 1:k){
    n=selected[,p]
    weakwei=weights[[p]][partialF[[p]]<F_threshold]
    if(IV_matrix[p,2]>1){
      weakwei=weightcal(weakwei)
      combineMajor[[p]]=cbind(G[n,selectmajor[[p]]],G[n,selectweak[[p]]]%*%weakwei)
    } 
    if(IV_matrix[p,2]<=1) combineMajor[[p]]=G[n,selectIV[[p]]]
  }
  #MAJOR
  
  hatX=rep(0,N)
  #Stage one: get estimated exposure
  for(p in 1:k){
    n=selected[,p]
    fit=lm(X[n]~combineMajor[[p]])
    hatX[n]=fit$fitted.values
  }
  return(list("hatX"=hatX,"selectmajor"=selectmajor,"selectIV"=selectIV))
}
#Assume no noise
MR_SPLIT_one_nonoise=function(X,Y,G,k=2,F_threshold=1000000){
  N=length(X)
  #split sample
  selected=NULL
  unselected=1:N
  J=ncol(G)
  n=sample(unselected,size=N/k,replace = FALSE)
  selected=cbind(selected,n)
  for(p in 1:k){
    if(p<k){
      n=sample(unselected[-selected],size=N/k,replace = FALSE)
      selected=cbind(selected,n)
    }
  }
  IV_matrix=matrix(0,ncol=4,nrow=2)
  
  selectIV=list()
  weights=list()
  
  for(p in 1:k){
    n=selected[,p]
    selectIV[[p]]=1:ncol(G)
    
    fit=lm(X[-n]~G[-n,])
    weights[[p]]=coef(fit)[-1]
    IV_matrix[p,3]=length(selectIV[[p]])
  }
  
  
  #select major and weak based on partial F statistics
  partialF=list()
  for(p in 1:k){
    n=selected[,p]
    full=lm(X[-n]~G[-n,selectIV[[p]]])
    Fs=NULL
    if(length(selectIV[[p]])==1){
      Fs=summary(lm(X[-n]~G[-n,selectIV[[p]]]))$fstat[1]
    }
    if(length(selectIV[[p]])>1){
      for(q in 1:length(selectIV[[p]])){
        reduced=lm(X[-n]~G[-n,selectIV[[p]][-q]])
        Fs[q]=anova(reduced,full)$F[2]
      }
    }
    partialF[[p]]=Fs
  }
  
  selectmajor=list()
  selectweak=list()
  for(p in 1:k){
    if(is.na(selectIV[p])==FALSE){
      selectmajor[[p]]=selectIV[[p]][partialF[[p]]>=F_threshold]
      selectweak[[p]]=selectIV[[p]][partialF[[p]]<F_threshold]
      IV_matrix[p,1]=sum(partialF[[p]]>=F_threshold)
      IV_matrix[p,2]=sum(partialF[[p]]<F_threshold)
    }
  }
  
  #Get coef
  combineMajor=list()
  for(p in 1:k){
    n=selected[,p]
    weakwei=weights[[p]][partialF[[p]]<F_threshold]
    if(IV_matrix[p,2]>1){
      weakwei=weightcal(weakwei)
      combineMajor[[p]]=cbind(G[n,selectmajor[[p]]],G[n,selectweak[[p]]]%*%weakwei)
    } 
    if(IV_matrix[p,2]<=1) combineMajor[[p]]=G[n,selectIV[[p]]]
  }
  #MAJOR
  
  hatX=rep(0,N)
  #Stage one: get estimated exposure
  for(p in 1:k){
    n=selected[,p]
    fit=lm(X[n]~combineMajor[[p]])
    hatX[n]=fit$fitted.values
  }
  return(list("hatX"=hatX,"selectmajor"=selectmajor,"selectIV"=selectIV))
}


#generate samples
generate_SNP_discrete=function(J,N,maf){
  Sigma = 0.8*diag(1,J,J)
  # 生成多变量正态分布数据 (N x J)
  continuous_data <- mvrnorm(n = N, mu = rep(0, J), Sigma = Sigma)
  
  # 初始化SNP矩阵
  SNP_data <- matrix(0, nrow = N, ncol = J)
  
  # 计算阈值用于将连续值转为0, 1, 2
  thresholds <- qnorm(c((1-maf)^2, 1-maf^2), sd=0.8^0.5,lower.tail = TRUE)
  
  # 根据固定的 MAF 将连续数据离散化为 0, 1, 2
  SNP_data <- ifelse(continuous_data < thresholds[1], 0, 
                     ifelse(continuous_data < thresholds[2], 1, 2))
  
  
  return(SNP_data)
}
generate_SNP=function(J,N){
  Sigma = 0.8*diag(1,J,J)
  SNP=mvrnorm(n=N,mu=rep(0,J),Sigma=Sigma)
  return(SNP)
}
generate_sample=function(beta,alpha_1,alpha_2,N,J,Sigma){
  J1=length(alpha_1)
  G=generate_SNP(J,N)
  error=mvrnorm(N,mu=c(0,0),Sigma=Sigma)
  #randomly choose SNPs that generates the esposure
  real=sample(c(1:J),size=J1,replace = FALSE)
  Greal=G[,real]
  #generate exposure
  X=as.matrix(Greal)%*%as.matrix(alpha_1)+error[,1]
  #generate outcome
  Y=X*beta+as.matrix(Greal)%*%as.matrix(alpha_2)+error[,2]
  return(list("X"=X,"Y"=Y,"G"=G,"real_IV"=real))
}

#other functions
identify_components_percentage <- function(A, valid_codes, invalid_codes) {
  # 创建一个向量，用于存储结果
  results <- character(length(A))
  
  # 循环遍历向量A
  for (i in seq_along(A)) {
    if (A[i] %in% valid_codes) {
      results[i] <- "valid"
    } else if (A[i] %in% invalid_codes) {
      results[i] <- "invalid"
    } else {
      results[i] <- "noise"
    }
  }
  
  # 将结果转换为因子，并设置水平
  results_factor <- factor(results, levels = c("valid", "invalid", "noise"))
  
  # 计算每个成分的数量
  counts <- table(results_factor)
  
  return(counts)
}
weightcal=function(weight){
  len=length(weight)
  a=(weight>=0)
  b=abs(weight)/sum(abs(weight))
  c=NULL
  for(i in 1:len){
    if(a[i])c[i]=b[i]
    if(!a[i])c[i]=-1*b[i]
  }
  return(c)
}

