library(readxl)
library(stringr)
library(rlang)
library(e1071)
library(glmnet)
library(pROC)

library(remotes); install_github("Mirrelijn/ecpc/Rpackage")
library(ecpc); ?ecpc
library(MASS)
library(penalized)
library(glmnet)
library(mvtnorm)
library(gglasso)
library(mgcv)
library(CVXR)
library(GRridge)
library(randomForest)
library(expm)
library(Rsolnp)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggraph)
library(igraph)
library(RColorBrewer)


remotes::install_github("DennisBeest/CoRF")
library(CoRF)
library(randomForestSRC)
library(scam)

#source('spatstat_vectra.R')




## FUNCTIONS FROM SPATSTAT_VECTRA

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





##### DATA PREPROCESSING ####

##### read feature matrix and clinical data
# clinical_raw <- read_excel("C:/Users/t.brug/Documents/Bosch, Erik/2019_06_20_Full spreadsheet HO105_cleaned.xls")
# features_raw = readRDS(path.expand('C:/Users/t.brug/Documents/Bosch, Erik/Extracted_Features_4thAug.RDS'))

clinical_raw <- read_excel("~/Studie/Thesis - Local/2019_06_20_Full spreadsheet HO105_cleaned.xls")
features_raw = readRDS(path.expand('~/Studie/Thesis - Local/Extracted_Features_4thAug.RDS'))


##### paste 'Centered_' in centered column name statistics for better regex finding in making groups
index_pcf = which(str_detect(colnames(features_raw),'(?<![Normalized_])pcf_radius'))
colnames(features_raw)[index_pcf] = paste0('Centered_', colnames(features_raw)[index_pcf])
index_other_stats = which(str_detect(colnames(features_raw),'(?<![Normalized_])[A-Z]_radius|(?<![Normalized_])[A-Z]dot_radius'))
colnames(features_raw)[index_other_stats] = paste0('Centered_', colnames(features_raw)[index_other_stats])

features_raw = features_raw[,c(colnames(features_raw)[1:86],sort(colnames(features_raw)[87:length(colnames(features_raw))]))]

##### make correct patient names and filter without Jstatistic
rownames(features_raw) <- str_replace_all(rownames(features_raw), pattern = '\\_[0-9]\\_', replacement = '\\_')
cols_without_J = str_subset(colnames(features_raw),'J_radius', negate = TRUE)
cols_without_Jdot = str_subset(colnames(features_raw),'Jdot_radius', negate = TRUE)
cols_without_J_Total = intersect(cols_without_J, cols_without_Jdot)

features_raw = features_raw[,cols_without_J_Total]


##### convert raw outcome on raw clinical data to binary outcome by thresholding
outcome_raw = rep(0,dim(clinical_raw)[1])
names(outcome_raw) = clinical_raw[['case_id']]

# set threshold months to 12, 9 or 6 months
threshold_months = 12


for (i in seq_along(outcome_raw)){
  if (clinical_raw[i, 'Event free survival [m]'] > threshold_months){
    outcome_raw[i] = 0
  } else if (clinical_raw[i,'EFS [0,1]'] == 1){
    outcome_raw[i] = 1
  } else {
    outcome_raw[i] = NA
  }
}

##### sync data with the patient MSI selected by biologist
image_check_path <- path.expand('~/Studie/Thesis - Local/spatstat_vectra/HO105 tumor_border annotation.xlsx')

HO105_image_check <- read_excel(image_check_path)


### Borders ###

# replace string 'NA' with NA in border column
HO105_image_check$`Border images`[HO105_image_check$`Border images` == 'NA'] = rep(NA)
# get the MSI's through regex finding: open bracket->MSI numbers->comma->MSI numbers-> closed bracket
border_images_unfil = str_extract_all(HO105_image_check$`Border images`,'\\[[0-9]+,[0-9]+\\]')
# create filtered border_images and HOnrs_border
border_images = border_images_unfil[!is.na(border_images_unfil)]
HOnrs_border = HO105_image_check$`HO105 nr`[!is.na(border_images_unfil)]
names(border_images) = HOnrs_border
# create string with sample names of tumor images that are checked by marit
border_images_checked = unlist(lapply(seq_along(HOnrs_border),function(x) paste0(HOnrs_border[x],'_',unlist(border_images[HOnrs_border][x]))))

warning('these border images are missing in features matrix ', setdiff(border_images_checked,rownames(features_raw)), ' removing from borders to evaluate')
if (!is_empty(setdiff(border_images_checked,rownames(features_raw)))){
  border_images_checked = intersect(border_images_checked,rownames(features_raw))
}

features_border_checked = features_raw[border_images_checked,]

##### averaging features for every  patient in the feature matrix with specific NA handling
features_border_checked = statisticPerPatient(features_border_checked, statistic = 'mean', na.handler = 'complete_cases')



##### find complete patient cases in outcome and features
patientID_outcome = names(outcome_raw[complete.cases(outcome_raw)])
patientID_features = rownames(features_border_checked[complete.cases(features_border_checked),])
patientID_complete = intersect(patientID_outcome, patientID_features)


##### select complete patient cases outcome and features in correct syntax for prediction
outcome_complete = outcome_raw[patientID_complete]
features_complete = features_border_checked[patientID_complete,]
patientID_complete = patientID_complete



##### FEATURE PREPROCESSING #####

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


#final features for prediction!
features = scale(deskew(delete_constant_features(features_complete), 1))
outcome = outcome_complete


###### RIDGE LOGISTIC REGRESSION ######

set.seed(0)
myalpha = 0 #alpha=0 ridge, alpha=1 LASSO
penaltyfactor = rep(1,ncol(features)) #1 is penalized, 0 is unpenalized
train_percentage = 0.8 
nrepeats = 250

aucs = numeric(nrepeats) 
indices0 = (1:length(outcome))[outcome==0]
indices1 = (1:length(outcome))[outcome==1]
for (k in 1:nrepeats) {
  ind0 = sample(indices0)
  ind1 = sample(indices1)
  indices_train = c(ind0[1:round(length(ind0)*train_percentage)],
                    ind1[1:round(length(ind1)*train_percentage)])
  indices_test = c(ind0[(round(length(ind0)*train_percentage)+1):length(ind0)],
                   ind1[(round(length(ind1)*train_percentage)+1):length(ind1)])
  best_lambda = cv.glmnet(features[indices_train,], outcome[indices_train], alpha=myalpha, penalty.factor = penaltyfactor, nfolds=5, family="binomial")$lambda.min
  model_train = glmnet(features[indices_train,], outcome[indices_train], alpha=myalpha, penalty.factor = penaltyfactor, family="binomial", lambda=best_lambda)
  prob_test = as.vector(predict(model_train, type="response", newx=features[indices_test,]))
  myroc = pROC::roc(outcome[indices_test]~prob_test,levels=c(0,1), direction="<", quiet=TRUE)
  aucs[k] = myroc$auc
  
}
{
  cat('threshold months is set to',threshold_months,fill = T)
  cat('mean of aucs',fill = T)
  cat(mean(aucs),fill = T)
  cat('sd of aucs',fill = T)
  cat(sd(aucs),fill = T)
  cat('2.5 and 97.5 quantiles of aucs',fill = T)
  cat(quantile(aucs,0.025),fill = T)
  cat(quantile(aucs,0.975),fill = T)#cross validated AUC
}
# threshold months is set to 6
# mean of aucs
# 0.5666111
# sd of aucs
# 0.1707911
# 2.5 and 97.5 quantiles of aucs
# 0.2284722
# 0.8611111


# threshold months is set to 12
# mean of aucs
# 0.5699429
# sd of aucs
# 0.1746525
# 2.5 and 97.5 quantiles of aucs
# 0.2571429
# 0.9364286


###### RIDGE LOGISTIC REGRESSION WITH GROUP STRUCTURE (ECPC) ######


# function for prediction depending on mygrouping
grouping_prediction <- function(outcome, features, mygrouping,threshold_months = threshold_months) {
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
  
  train_percentage = 0.8 
  nrepeats = 50 
  
  aucs = numeric(nrepeats)
  gammas = matrix(nrow=nrepeats, ncol=ngroups) #group weights
  indices0 = (1:length(outcome))[outcome==0]
  indices1 = (1:length(outcome))[outcome==1]
  
  for (k in 1:nrepeats) { 
    print(paste("REPEAT", k))
    ind0 = sample(indices0)
    ind1 = sample(indices1)
    indices_train = c(ind0[1:round(length(ind0)*train_percentage)],
                      ind1[1:round(length(ind1)*train_percentage)])
    indices_test = c(ind0[(round(length(ind0)*train_percentage)+1):length(ind0)],
                     ind1[(round(length(ind1)*train_percentage)+1):length(ind1)])
    fit = ecpc(Y=outcome[indices_train], X=features[indices_train,], Y2=outcome[indices_test], X2=features[indices_test,], 
               groupings=mygrouping, model="logistic", postselection=FALSE)
    myroc = pROC::roc(outcome[indices_test]~fit$Ypred,levels=c(0,1), direction="<", quiet=TRUE)
    aucs[k] = myroc$auc
    gammas[k,] = fit$gamma
  }
  cat('threshold months is set to',threshold_months,fill = T)
  cat('mean of aucs',fill = T)
  cat(mean(aucs),fill = T)
  cat('sd of aucs',fill = T)
  cat(sd(aucs),fill = T)
  cat('2.5 and 97.5 quantiles of aucs',fill = T)
  cat(quantile(aucs,0.025),fill = T)
  cat(quantile(aucs,0.975),fill = T) #cross validated AUC
  cat('estimated group weights',fill = T)
  cat(colMeans(gammas)/sum(colMeans(gammas)),fill = T) #estimated group weights
  
  return(list('aucs' = aucs,'gammas' = gammas,'fit' = fit, 'myroc' = myroc))
}

##### groups #####
## 0
GRallfeatgroup = 1:1366
## 1
GRsimple_spatial <- 1:24
GRcomplex_spatial <- 25:1366
# 2
GRcounts = c(1:32) # L32 Regex: which(str_detect(colnames(features),'counts_[a-z]+|density|X2stat')) 
GRmed_mad_ss = c(33:86) # L54 Regex: which(str_detect(colnames(features),'MED|MAD|distance_ratio'))
GRcentered = c(87:726) # L640 Regex: which(str_detect(colnames(features), 'Centered_'))
GRnormalized = c(727:1366) # L640 Regex: which(str_detect(colnames(features), 'Normalized_'))
## 3
#GRcounts
GRmed_ss = c(33:58,85,86) # L28 Regex: which(str_detect(colnames(features),'MED|distance_ratio'))
GRmad = c(59:84) # L28 Regex: which(str_detect(colnames(features),'MAD'))
grFstat = c(87:118,727:758) # L64 which(str_detect(colnames(features),'F_radius'))
grGstat = c(119:278,759:918) # L320 which(str_detect(colnames(features),'G_radius|Gdot_radius'))
grKstat = c(279:438, 919:1078)  # L320 which(str_detect(colnames(features),'K_radius|Kdot_radius'))
grLstat = c(439:598, 1079:1238) # L320 which(str_detect(colnames(features),'L_radius|Ldot_radius'))
GRpcfstat = c(599:726, 1239:1366) # L256 which(str_detect(colnames(features),'pcf_radius'))
## 4
#GRcounts 
#GRmed_ss
#GRmad
grFstat_close = which(str_detect(colnames(features),'F_radius [0-9]_')) # L24
grFstat_far = which(str_detect(colnames(features),'F_radius [0-9][0-9]_')) # L40 
grGstat_close = which(str_detect(colnames(features),'G_radius [0-9]_|Gdot_radius [0-9]_')) # L120 
grGstat_far = which(str_detect(colnames(features),'G_radius [0-9][0-9]_|Gdot_radius [0-9][0-9]_')) # L200 
grKstat_close = which(str_detect(colnames(features),'K_radius [0-9]_|Kdot_radius [0-9]_')) # L120 
grKstat_far = which(str_detect(colnames(features),'K_radius [0-9][0-9]_|Kdot_radius [0-9][0-9]_')) # L200 
grLstat_close = which(str_detect(colnames(features),'L_radius [0-9]_|Ldot_radius [0-9]_')) # L120 
grLstat_far = which(str_detect(colnames(features),'L_radius [0-9][0-9]_|Ldot_radius [0-9][0-9]_')) # L200 
GRpcfstat_close = which(str_detect(colnames(features),'pcf_radius [0-9]_')) # L96 
GRpcfstat_far = which(str_detect(colnames(features),'pcf_radius [0-9][0-9]_')) # L160
# 5
grMacrophage = which(str_detect(colnames(features),'Macrophage|distance_ratio')) #L547
grTcells = which(str_detect(colnames(features),'Tcells|distance_ratio')) #L547
grTumors = which(str_detect(colnames(features),'Tumor|distance_ratio')) #L547
grOthers = which(str_detect(colnames(features),'Others')) #L545

##### mygroupings #####
mygrouping0 = list(GRallfeatgroup)
mygrouping1 = list(GRsimple_spatial,GRcomplex_spatial)
mygrouping2 = list(GRcounts,GRmed_mad_ss,GRcentered,GRnormalized) # L1366 sum(unlist(lapply(mygrouping,function(x) length(x))))
mygrouping3 = list(GRcounts,GRmed_ss,GRmad,grFstat,grGstat,grKstat,grLstat,GRpcfstat) # L1366 sum(unlist(lapply(mygrouping,function(x) length(x))))
mygrouping4 = list(GRcounts,GRmed_ss,GRmad,
                   grFstat_close,grFstat_far,grGstat_close, grGstat_far,
                   grKstat_close,grKstat_far,grLstat_close,grLstat_far,
                   GRpcfstat_close,GRpcfstat_far) # L1366 sum(unlist(lapply(mygrouping,function(x) length(x))))
mygrouping5 = list(grMacrophage,grTcells,grTumors,grOthers) # L1962 sum(unlist(lapply(mygrouping,function(x) length(x))))
mygrouping35 = list(mygrouping3,mygrouping5)

##### RIDGE LOGISTIC REGRESSION WITH GROUP STRUCTURE (ECPC) 6 months #####
set.seed(0)
# ridge_outg0m6 = grouping_prediction(outcome, features, mygrouping0,threshold_months) # make sure to change threshold_months
# ridge_outg1m6 = grouping_prediction(outcome, features, mygrouping1,threshold_months)
# ridge_outg2m6 = grouping_prediction(outcome, features, mygrouping2,threshold_months)
ridge_outg3m6 = grouping_prediction(outcome, features, mygrouping3,threshold_months)
saveRDS(object = ridge_outg3m6,file = 'ridge_outg3m6.RDS')
# [1] "REPEAT 8"
# [1] "Estimate global tau^2 (equiv. global ridge penalty lambda)"
# Error in qr.coef(qr(a, LAPACK = TRUE), b) : 
#   error code 42 from Lapack routine 'dtrtrs'
# In addition: Warning message:
#   In ecpc(Y = outcome[indices_train], X = features[indices_train,  :
#                                                      Cross-validated global penalty lambda was >10^12 and set to 100
#                                                    Error: Matrix inversion failed. Please increase lambda1 and/or lambda2 
# ridge_outg4m6 = grouping_prediction(outcome, features, mygrouping4,threshold_months)
ridge_outg5m6 = grouping_prediction(outcome, features, mygrouping5,threshold_months)
saveRDS(object = ridge_outg5m6,file = 'ridge_outg5m6.RDS')
# threshold months is set to 6
# mean of aucs
# 0.6533333
# sd of aucs
# 0.1250145
# 2.5 and 97.5 quantiles of aucs
# 0.4722222
# 0.9104167
# estimated group weights
# 0.2517929 0.006555848 0.370991 0.3706602
# Warning message:
#   In ecpc(Y = outcome[indices_train], X = features[indices_train,  :
#                                                      Cross-validated global penalty lambda was >10^12 and set to 100
ridge_outg35m6 = grouping_prediction(outcome, features, mygrouping35,threshold_months)
saveRDS(object = ridge_outg35m6,file = 'ridge_outg35m6.RDS')
# [1] "REPEAT 9"
# [1] "Estimate global tau^2 (equiv. global ridge penalty lambda)"
# Error in qr.coef(qr(a, LAPACK = TRUE), b) : 
#   error code 39 from Lapack routine 'dtrtrs'
# Error: Matrix inversion failed. Please increase lambda1 and/or lambda2 


stop('rerun script to continue with new outcome with new threshold month and therefore outcome')

##### RIDGE LOGISTIC REGRESSION WITH GROUP STRUCTURE (ECPC) 12 months #####
set.seed(0)
# ridge_outg0m12 = grouping_prediction(outcome, features, mygrouping0,threshold_months)
# ridge_outg1m12 = grouping_prediction(outcome, features, mygrouping1,threshold_months)
# ridge_outg2m12 = grouping_prediction(outcome, features, mygrouping2,threshold_months)
ridge_outg3m12 = grouping_prediction(outcome, features, mygrouping3,threshold_months)
saveRDS(object = ridge_outg3m12,file = 'ridge_outg3m12.RDS')
# threshold months is set to 12
# mean of aucs
# 0.6714286
# sd of aucs
# 0.1540512
# 2.5 and 97.5 quantiles of aucs
# 0.4064286
# 0.9364286
# estimated group weights
# 0 0.7482927 0.0340254 0.1426701 0.002286172 0.001657105 0.04497709 0.02609149
# ridge_outg4m12 = grouping_prediction(outcome, features, mygrouping3,threshold_months)
ridge_outg5m12 = grouping_prediction(outcome, features, mygrouping5,threshold_months)
saveRDS(object = ridge_outg5m12,file = 'ridge_outg5m12.RDS')
# threshold months is set to 12
# mean of aucs
# 0.6434286
# sd of aucs
# 0.1609261
# 2.5 and 97.5 quantiles of aucs
# 0.2921429
# 0.9235714
# estimated group weights
# 0.3426152 0.004473166 0.2876232 0.3652884
# Warning messages:
#   1: In ecpc(Y = outcome[indices_train], X = features[indices_train,  :
#                                                         Cross-validated global penalty lambda was >10^12 and set to 100
#   2: In ecpc(Y = outcome[indices_train], X = features[indices_train,  :
#                                                         Cross-validated global penalty lambda was >10^12 and set to 100
#   3: In ecpc(Y = outcome[indices_train], X = features[indices_train,  :
#                                                         Cross-validated global penalty lambda was >10^12 and set to 100
ridge_outg35m12 = grouping_prediction(outcome, features, mygrouping35,threshold_months)
saveRDS(object = ridge_outg35m12,file = 'ridge_outg35m12.RDS')
# threshold months is set to 12
# mean of aucs
# 0.6297143
# sd of aucs
# 0.1775355
# 2.5 and 97.5 quantiles of aucs
# 0.2764286
# 0.8857143
# estimated group weights
# 8.332623e-05 0.6668872 0.02658449 0.1603943 0.008039079 0.001159193 0.03955409 0.03937858 0.0234668 0.0004502759 0.01358855 0.02041418


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

##### grouping 3 for Random Forest ####### 
group3_per_feature = numeric(1366)
group3_per_feature[GRcounts] = 1
group3_per_feature[GRmed_ss] = 2
group3_per_feature[GRmad] = 3
group3_per_feature[grFstat] = 4
group3_per_feature[grGstat] = 5
group3_per_feature[grKstat] = 6
group3_per_feature[grLstat] = 7
group3_per_feature[GRpcfstat] = 8
##### grouping 5 for Random Forest ####### 
group5_per_feature = numeric(1366)
group5_per_feature[grMacrophage] = 1
group5_per_feature[grTcells] = 2
group5_per_feature[grTumors] = 3
group5_per_feature[grOthers] = 4



##### RANDOM FOREST 6 months #####
# threshold months is set to 6
print(paste('threshold months is set to',threshold_months))
set.seed(0)
DF <- data.frame(outcome, features)
Forest <- rfsrc(outcome ~ .,data=DF,ntree=20000,var.used="all.trees",importance="TRUE",nodesize=2,seed=1)

roc_Forest = pROC::roc(outcome~Forest$predicted.oob, levels=c(0,1), direction="<")
plot(roc_Forest)
roc_Forest$auc #Out Of Bag performance
# Area under the curve: 0.6398
# Forest$importance #feature importances
rf_outm6 = list('Forest' = Forest, 'roc_Forest' = roc_Forest)
saveRDS(object = rf_outm6,file = 'rf_outm6.RDS')

##### RANDOM FOREST WITH GROUP STRUCTURE (CORF) 6 months  ####### 

set.seed(0)
rf_outg3m6 = reforest(outcome, features, Forest, group3_per_feature, threshold_months)
# threshold months is set to 6
# Area under curve:  0.6505376

# save outcome
saveRDS(object = rf_outg3m6,file = 'rf_outg3m6.RDS')

rf_outg5m6 = reforest(outcome, features, Forest, group5_per_feature, threshold_months)
# threshold months is set to 6
# Area under curve:  0.6408602

# save outcome
saveRDS(object = rf_outg5m6,file = 'rf_outg5m6.RDS')



stop('rerun script to continue with new outcome with new threshold month and therefore outcome')

##### RANDOM FOREST 12 months #####
# threshold months is set to 12
print(paste('threshold months is set to',threshold_months))
set.seed(0)
DF <- data.frame(outcome, features)
Forest <- rfsrc(outcome ~ .,data=DF,ntree=20000,var.used="all.trees",importance="TRUE",nodesize=2,seed=1)

roc_Forest = pROC::roc(outcome~Forest$predicted.oob, levels=c(0,1), direction="<")
plot(roc_Forest)
roc_Forest$auc #Out Of Bag performance
# Area under the curve: 0.6917
# Forest$importance #feature importances
rf_outm12 = list('Forest' = Forest, 'roc_Forest' = roc_Forest)
saveRDS(object = rf_outm12,file = 'rf_outm12.RDS')

##### RANDOM FOREST WITH GROUP STRUCTURE (CORF) 12 months  ####### 

set.seed(0)
rf_outg3m12 = reforest(outcome, features, Forest, group3_per_feature, threshold_months)
# threshold months is set to 12
# Area under curve:  0.6917211

# save outcome
saveRDS(object = rf_outg3m12,file = 'rf_outg3m12.RDS')

rf_outg5m12 = reforest(outcome, features, Forest, group5_per_feature, threshold_months)
# threshold months is set to 12
# Area under curve:  0.6753813

# save outcome
saveRDS(object = rf_outg5m12,file = 'rf_outg5m12.RDS')
