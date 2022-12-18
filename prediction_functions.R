
library(rlang)
library(e1071)
library(glmnet)
library(pROC)

library(remotes); install_github("Mirrelijn/ecpc/Rpackage")
library(ecpc); #?ecpc
library(MASS)
library(penalized)
library(mvtnorm)
library(gglasso)
library(mgcv)
#library(CVXR)
library(GRridge)
library(randomForest)
library(expm)
library(Rsolnp)
library(dplyr)
library(purrr)
library(foreach)
library(doMC)

#remotes::install_github("DennisBeest/CoRF")
library(CoRF)
library(randomForestSRC)
library(scam)



statisticPerPatient <- function(mat, statistic = 'mean', na.handler = 'complete_cases'){
  
  # browser()
  if (isTRUE(na.handler == 'complete_cases')){
    # cat('samplenames with missing values removed for analyse:', rownames(mat[!complete.cases(mat),]), '. if empty either 0 or 1 samplenames are removed')
    na.rm = FALSE # default setting
    mat = mat[complete.cases(mat),] # complete cases
  } else if (isTRUE(na.handler == 'ignore_na')){
    na.rm = TRUE
  } else {
    stop('input na.handler must be either complete_cases or ignore_na')
  }
  
  if (isFALSE(statistic %in% c('mean','median'))){
    stop('input statistic must be either mean or median')
  }
  
  samplenames = rownames(mat)
  prediction_statistics = colnames(mat)
  nrs = unique(str_remove_all(samplenames,pattern = '\\_\\[[0-9]+,[0-9]+\\]'))
  
  mat_allpatients = matrix(NA, nrow = length(nrs), ncol = length(prediction_statistics))
  rownames(mat_allpatients) = nrs
  colnames(mat_allpatients) = prediction_statistics
  # browser()
  for (patientnr in nrs){
    
    samplenames_patient = str_subset(samplenames, patientnr)
    if (isTRUE(length(samplenames_patient) == 1) ){
      mat_patient = t(as.matrix(mat[samplenames_patient,]))
    } else {
      mat_patient = mat[samplenames_patient,]
    }
    mat_allpatients[patientnr,] = apply(mat_patient, 2, statistic, na.rm = na.rm)
  }
  mat_allpatients[is.nan(mat_allpatients)] = NA
  
  return(mat_allpatients)
}



#delete constant features
delete_constant_features <- function(myfeatures) {
  constant_features = c()
  for (j in 1:ncol(myfeatures)) {
    if (max(myfeatures[,j]) - min(myfeatures[,j]) == 0) {
      constant_features = c(constant_features,j)
    }
  }
  if (is.null(constant_features)) {
    return(myfeatures)
  } else {
    return(myfeatures[,-constant_features])
  }
}


#log transform of skewed features
deskew <- function(myfeatures, threshold) {
  myfeatures_deskewed = myfeatures
  for (j in 1:ncol(myfeatures)) {
    if (skewness(myfeatures[,j])>threshold & sum(myfeatures[,j]<=0)==0) { #skewed and positive
      for (i in 1:nrow(myfeatures)) {
        myfeatures_deskewed[i,j] = log(myfeatures[i,j])
      }
    }
  }
  return(myfeatures_deskewed)
}


# function for prediction depending on mygrouping
grouping_prediction <- function(outcome, features, mygrouping,
                                nrepeats = 250, train_percentage = 0.8,
                                ncore = 10, maxsel=100) {
  registerDoMC(ncore)
  set.seed(0)
  ndepth = vec_depth(mygrouping)
  if (ndepth > 2){
    # multiple groupings
    ngroupings = length(mygrouping)
    ngroups = 0
    for (groupingindex in 1:ngroupings){
      ngroups = ngroups + length(mygrouping[[groupingindex]])
    }
  } else{
    # one grouping
    ngroupings = 1
    ngroups = length(mygrouping)
    mygrouping = list(mygrouping)
  }
 
  aucs = numeric(nrepeats)
  gammas = matrix(nrow=nrepeats, ncol=ngroups) #group weights
  indices0 = (1:length(outcome))[outcome==0]
  indices1 = (1:length(outcome))[outcome==1]
  
  out = foreach(k=1:nrepeats) %dopar% {
    ind0 = sample(indices0)
    ind1 = sample(indices1)
    indices_train = c(ind0[1:round(length(ind0)*train_percentage)],
                      ind1[1:round(length(ind1)*train_percentage)])
    indices_test = c(ind0[(round(length(ind0)*train_percentage)+1):length(ind0)],
                     ind1[(round(length(ind1)*train_percentage)+1):length(ind1)])
    fit = ecpc(Y=outcome[indices_train], X=features[indices_train,], Y2=outcome[indices_test], X2=features[indices_test,], 
               groupings=mygrouping, model="logistic", maxsel=maxsel)
    myroc = pROC::roc(outcome[indices_test]~fit$Ypred,levels=c(0,1), direction="<", quiet=TRUE)
    return(list(auc=myroc$auc, gamma=fit$gamma, fit=fit))
  }
  fits = list()
  for (k in 1:nrepeats) { 
    aucs[k] = out[[k]]$auc
    gammas[k,] = out[[k]]$gamma
    fits[[k]] = out[[k]]$fit
  }
  cat('mean of aucs',fill = T)
  cat(mean(aucs),fill = T)
  cat('sd of aucs',fill = T)
  cat(sd(aucs),fill = T)
  cat('2.5 and 97.5 quantiles of aucs',fill = T)
  cat(quantile(aucs,0.025),fill = T)
  cat(quantile(aucs,0.975),fill = T) #cross validated AUC
  cat('estimated group weights',fill = T)
  cat(colMeans(gammas)/sum(colMeans(gammas)),fill = T) #estimated group weights

  finalfit = ecpc(Y=outcome, X=features,
                  Y2=outcome, X2=features,
               groupings=mygrouping, model="logistic", maxsel=maxsel)


  return(list('aucs' = aucs,'gammas' = gammas, 'fits' = fits, finalfit=finalfit))
}


##### RANDOM FOREST ######
reforest <- function(outcome, features, Forest, group_per_feature, threshold_months){
  set.seed(0)
  DF <- data.frame(outcome, features)
  
  #number of times each feature is used in the ordinary random forest
  VarUsed = Forest$var.used
  
  nfeatures = ncol(features)
  preds2 = c()
  for (j in 1:nfeatures) {
    g = group_per_feature[j]
    pred = mean(VarUsed[group_per_feature==g]) / sum(VarUsed)
    preds2[j] = max(pred-1/nfeatures,0)
  }
  
  #number of features randomly selected as candidates for splitting a node
  Mtry <- ceiling(sqrt(sum(preds2!=0))) 
  
  #run CORF random forest
  RefittedCoRF <- rfsrc(outcome ~ .,data=DF,ntree=20000,var.used="all.trees",importance="TRUE",
                        xvar.wt=preds2,mtry=Mtry,nodesize=2,setseed=1)
  
  #Out Of Bag performance
  roc_corf = pROC::roc(outcome~RefittedCoRF$predicted.oob, levels=c(0,1), direction="<")
  #plot(roc_corf)
  # feature importances
  #RefittedCoRF$importance 
  
  cat('threshold months is set to',threshold_months,fill = T)
  cat('Area under curve: ',roc_corf$auc,fill = T)
  
  return(list('RefittedCoRF' = RefittedCoRF, 'roc_corf' = roc_corf))
}


