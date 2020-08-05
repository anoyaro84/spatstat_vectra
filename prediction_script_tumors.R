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
library(ggplot2)
library(ggraph)
library(igraph)
library(RColorBrewer)

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





## DATA PREPROCESSING

# step 1 read feature matrix and clinical data
clinical_raw <- read_excel("C:/Users/t.brug/Documents/Bosch, Erik/2019_06_20_Full spreadsheet HO105_cleaned.xls")
features_raw = readRDS(path.expand('C:/Users/t.brug/Documents/Bosch, Erik/Extracted_Features_4thAug.RDS'))

clinical_raw <- read_excel("~/Studie/Thesis - Local/2019_06_20_Full spreadsheet HO105_cleaned.xls")
features_raw = readRDS(path.expand('~/Studie/Thesis - Local/Extracted_Features_4thAug.RDS'))


# paste 'Centered_' in centered column name statistics for better regex finding in making groups
index_pcf = which(str_detect(colnames(features_raw),'(?<![Normalized_])pcf_radius'))
colnames(features_raw)[index_pcf] = paste0('Centered_', colnames(features_raw)[index_pcf])
index_other_stats = which(str_detect(colnames(features_raw),'(?<![Normalized_])[A-Z]_radius|(?<![Normalized_])[A-Z]dot_radius'))
colnames(features_raw)[index_other_stats] = paste0('Centered_', colnames(features_raw)[index_other_stats])

features_raw = features_raw[,c(colnames(features_raw)[1:86],sort(colnames(features_raw)[87:length(colnames(features_raw))]))]

# make correct patient names and filter without Jstatistic
rownames(features_raw) <- str_replace_all(rownames(features_raw), pattern = '\\_[0-9]\\_', replacement = '\\_')
cols_without_J = str_subset(colnames(features_raw),'J_radius', negate = TRUE)
cols_without_Jdot = str_subset(colnames(features_raw),'Jdot_radius', negate = TRUE)
cols_without_J_Total = intersect(cols_without_J, cols_without_Jdot)

features_raw = features_raw[,cols_without_J_Total]


# step 2 compute raw outcome on raw clinical data
outcome_raw = rep(0,dim(clinical_raw)[1])
names(outcome_raw) = clinical_raw[['case_id']]

# set threshold months to 12, 9 or 6 months
threshold_months = 9

for (i in seq_along(outcome_raw)){
  if (clinical_raw[i, 'Event free survival [m]'] > threshold_months){
    outcome_raw[i] = 0
  } else if (clinical_raw[i,'EFS [0,1]'] == 1){
    outcome_raw[i] = 1
  } else {
    outcome_raw[i] = NA
  }
}

image_check_path <- path.expand('~/Studie/Thesis - Local/spatstat_vectra/HO105 tumor_border annotation.xlsx')

HO105_image_check <- read_excel(image_check_path)
# get the MSI's through regex finding: open bracket->MSI numbers->comma->MSI numbers-> closed bracket
tumor_images_unfil = str_extract_all(HO105_image_check$`Tumor images`,'\\[[0-9]+,[0-9]+\\]')

# create filtered tumor_images and HOnrs
tumor_images = tumor_images_unfil[!sapply(tumor_images_unfil,is_empty)]
HOnrs = HO105_image_check$`HO105 nr`[!sapply(tumor_images_unfil,is_empty)]
names(tumor_images) = HOnrs
# create string with sample names of tumor images that are checked by marit
tumor_images_checked = unlist(lapply(seq_along(HOnrs),function(x) paste0(names(tumor_images)[x],'_',unlist(tumor_images[HOnrs][x]))))


warning('these are missing in features matrix ', setdiff(tumor_images_checked,rownames(features_raw)))
features_tumor_checked = features_raw[tumor_images_checked,]


# # step 3a Mean and Complete cases and with threshold
features_tumor_checked = statisticPerPatient(features_tumor_checked, statistic = 'mean', na.handler = 'complete_cases')


# step 4a filter on outcome NA's
patientID_outcome = names(outcome_raw[complete.cases(outcome_raw)])
patientID_features = rownames(features_tumor_checked[complete.cases(features_tumor_checked),])
patientID_complete = intersect(patientID_outcome, patientID_features)

# step 5a syntax input for prediction
outcome_complete = outcome_raw[patientID_complete]
patientID_complete = patientID_complete
features_complete = features_tumor_checked[patientID_complete,]
sum(clinical_raw$`Multiplex general panel` == 'yes', na.rm = TRUE) #89 patients multiplex
length(patientID_complete) # 86 patients viable for prediction


## FEATURE PREPROCESSING

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


## UNIVARIATE LOGISTIC REGRESSION

#univariate significant features

set.seed(0)
output = data.frame()
i = 1
for (j in 1:ncol(features)) {
  pvalue = summary(glm(outcome~features[,j], family = "binomial"))$coefficients[2,4]
  if (pvalue<0.001) {
    output[i,1] = j
    output[i,2] = colnames(features)[j]
    output[i,3] = pvalue
    i = i+1
  }
}
# View(output)
# saveRDS(object = output, file = 'prediction_tumors_uni.RDS')




## RIDGE LOGISTIC REGRESSION

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
mean(aucs); sd(aucs); quantile(aucs,0.025); quantile(aucs,0.975) #cross validated AUC
# [1] 0.6887429
# [1] 0.1145463
# 2.5% 
# 0.4714286 
# 97.5% 
# 0.8857143 




## RIDGE LOGISTIC REGRESSION WITH GROUP STRUCTURE (ECPC)


# function for prediction depending on mygrouping
grouping_prediction <- function(mygrouping) {
  set.seed(0)
  ngroups = length(mygrouping)
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
               groupings=list(mygrouping), model="logistic", postselection=FALSE)
    myroc = pROC::roc(outcome[indices_test]~fit$Ypred,levels=c(0,1), direction="<", quiet=TRUE)
    aucs[k] = myroc$auc
    gammas[k,] = fit$gamma
    
  }
  mean(aucs); sd(aucs); quantile(aucs,0.025); quantile(aucs,0.975) #cross validated AUC
  colMeans(gammas)/sum(colMeans(gammas)) #estimated group weights
  
  return(list(aucs,gammas))
}

# groupings per simple spatial (until density) vs spatia
mygrouping1 = list(1:24, 25:1366) 

out = grouping_prediction(mygrouping1)
aucs = out[[1]]
gammas = out[[2]]
mean(aucs); sd(aucs); quantile(aucs,0.025); quantile(aucs,0.975); #cross validated AUC
# [1] 0.728
# [1] 0.1170212
# 2.5% 
# 0.4667857 
# 97.5% 
# 0.8967857
colMeans(gammas)/sum(colMeans(gammas)) #estimated group weights
# [1] 0.93975114 0.06024886

# groupings per groups of statistics general (distinct groups)
counts = c(1:32) # L32 Regex: which(str_detect(colnames(features),'counts_[a-z]+|density|X2stat')) 
MED_MAD_SS = c(33:86) # L54 Regex: which(str_detect(colnames(features),'MED|MAD|distance_ratio'))
centered = c(87:726) # L640 Regex: which(str_detect(colnames(features), 'Centered_'))
normalized = c(727:1366) # L640 Regex: which(str_detect(colnames(features), 'Normalized_'))

mygrouping2 = list(counts,MED_MAD_SS,normalized,centered) # L1366 sum(unlist(lapply(mygrouping,function(x) length(x))))
# all features used: setdiff(colnames(features),unique(colnames(features)[unlist(mygrouping)])) = empty

out = grouping_prediction(mygrouping2)
aucs = out[[1]]
gammas = out[[2]]
mean(aucs); sd(aucs); quantile(aucs,0.025); quantile(aucs,0.975); #cross validated AUC
# [1] 0.7594286
# [1] 0.09643197
# 2.5% 
# 0.5778571 
# 97.5% 
# 0.9364286 
colMeans(gammas)/sum(colMeans(gammas)) #estimated group weights
# [1] 0.59404454 0.31990528 0.05679378 0.02925640


# groupings per statistics detailed (distinct groups)
counts_ss = c(1:32,85,86) # L34 Regex: which(str_detect(colnames(features),'counts_[a-z]+|density|X2stat|distance_ratio')) 
MED = c(33:58) # L26 Regex: which(str_detect(colnames(features),'MED'))
# to view: View(colnames(features)[which(str_detect(colnames(features),'MED'))])
MAD = c(59:84) # L26 Regex: which(str_detect(colnames(features),'MAD'))
Fstat = c(87:118,727:758) # L64 which(str_detect(colnames(features),'F_radius'))
Gstat = c(119:278,759:918) # L320 which(str_detect(colnames(features),'G_radius|Gdot_radius'))
Kstat = c(279:438, 919:1078)  # L320 which(str_detect(colnames(features),'K_radius|Kdot_radius'))
Lstat = c(439:598, 1079:1238) # L320 which(str_detect(colnames(features),'L_radius|Ldot_radius'))
pcfstat = c(599:726, 1239:1366) # L256 which(str_detect(colnames(features),'pcf_radius'))

mygrouping3 = list(counts_ss,MED,MAD,Fstat,Gstat,Kstat,Lstat,pcfstat) # L1366 sum(unlist(lapply(mygrouping,function(x) length(x))))
# all features used: setdiff(colnames(features),unique(colnames(features)[unlist(mygrouping)])) = empty

out = grouping_prediction(mygrouping3)
aucs = out[[1]]
gammas = out[[2]]
mean(aucs); sd(aucs); quantile(aucs,0.025); quantile(aucs,0.975); #cross validated AUC
# [1] 0.7642857
# [1] 0.1175641
# 2.5% 
# 0.5096429 
# 97.5% 
# 0.9428571 
colMeans(gammas)/sum(colMeans(gammas)) #estimated group weights
# [1] 0.21195463 0.53207727 0.05803904 0.10556036 0.02697113 0.01887658 0.01803074 0.02849024


# groupings per statistics detailed with radius close (<10) and far (>=10) (distinct groups)
counts_ss = c(1:32,85,86) # L34 Regex: which(str_detect(colnames(features),'counts_[a-z]+|density|X2stat|distance_ratio')) 
MED = c(33:58) # L26 Regex: which(str_detect(colnames(features),'MED'))
MAD = c(59:84) # L26 Regex: which(str_detect(colnames(features),'MAD'))
Fstat_close = which(str_detect(colnames(features),'F_radius [0-9]_')) # L24
Fstat_far = which(str_detect(colnames(features),'F_radius [0-9][0-9]_')) # L40 
Gstat_close = which(str_detect(colnames(features),'G_radius [0-9]_|Gdot_radius [0-9]_')) # L120 
Gstat_far = which(str_detect(colnames(features),'G_radius [0-9][0-9]_|Gdot_radius [0-9][0-9]_')) # L200 
Kstat_close = which(str_detect(colnames(features),'K_radius [0-9]_|Kdot_radius [0-9]_')) # L120 
Kstat_far = which(str_detect(colnames(features),'K_radius [0-9][0-9]_|Kdot_radius [0-9][0-9]_')) # L200 
Lstat_close = which(str_detect(colnames(features),'L_radius [0-9]_|Ldot_radius [0-9]_')) # L120 
Lstat_far = which(str_detect(colnames(features),'L_radius [0-9][0-9]_|Ldot_radius [0-9][0-9]_')) # L200 
pcfstat_close = which(str_detect(colnames(features),'pcf_radius [0-9]_')) # L96 
pcfstat_far = which(str_detect(colnames(features),'pcf_radius [0-9][0-9]_')) # L160 

mygrouping4 = list(counts_ss,MED,MAD,
                  Fstat_close,Fstat_far,Gstat_close, Gstat_far,
                  Kstat_close,Kstat_far,Lstat_close,Lstat_far,
                  pcfstat_close,pcfstat_far) # L1366 sum(unlist(lapply(mygrouping,function(x) length(x))))
# all features used: setdiff(colnames(features),unique(colnames(features)[unlist(mygrouping)])) = empty

out = grouping_prediction(mygrouping4)
aucs = out[[1]]
gammas = out[[2]]
mean(aucs); sd(aucs); quantile(aucs,0.025); quantile(aucs,0.975); #cross validated AUC
# [1] 0.7265714
# [1] 0.1047402
# 2.5% 
# 0.5175 
# 97.5% 
# 0.8825
colMeans(gammas)/sum(colMeans(gammas)) #estimated group weights
# [1] 2.595697e-01 2.447669e-01 8.870406e-02 7.155332e-03 1.665722e-01 5.757167e-06 9.541964e-02 7.450179e-02
# [9] 9.916802e-03 1.574040e-02 1.577482e-02 1.108556e-02 1.078710e-02


# groupings per phenotype (overlapping groups)
Macrophage = which(str_detect(colnames(features),'Macrophage|distance_ratio')) #L547
Tcells = which(str_detect(colnames(features),'Tcells|distance_ratio')) #L547
Tumors = which(str_detect(colnames(features),'Tumor|distance_ratio')) #L547
Others = which(str_detect(colnames(features),'Others')) #L545

mygrouping5 = list(Macrophage,Tcells,Tumors,Others) # L1962 sum(unlist(lapply(mygrouping,function(x) length(x))))
# all features used: setdiff(colnames(features),unique(colnames(features)[unlist(mygrouping)])) = empty

out = grouping_prediction(mygrouping5)
aucs = out[[1]]
gammas = out[[2]]
mean(aucs); sd(aucs); quantile(aucs,0.025); quantile(aucs,0.975); #cross validated AUC
# [1] 0.6817143
# [1] 0.116739
# 2.5% 
# 0.5142857 
# 97.5% 
# 0.9157143 
colMeans(gammas)/sum(colMeans(gammas)) #estimated group weights
# [1] 0.53283604 0.06991251 0.13806352 0.25918793



## RANDOM FOREST WITH GROUP STRUCTURE (CORF)

