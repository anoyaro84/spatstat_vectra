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

# take mean per patient for Tumors and change colnames accordingly
features_tumor_checked = statisticPerPatient(features_tumor_checked, statistic = 'mean', na.handler = 'complete_cases')
colnames(features_tumor_checked) <- paste0('Tumors_im_',colnames(features_tumor_checked))


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

# take mean per patient for Borders and change colnames accordingly
features_border_checked = statisticPerPatient(features_border_checked, statistic = 'mean', na.handler = 'complete_cases')
colnames(features_border_checked) <- paste0('Borders_im_',colnames(features_border_checked))



# merge features Tumors with features Borders
features_checked_df = merge(features_tumor_checked, features_border_checked, by = 'row.names', all = TRUE)
features_checked = as.matrix(features_checked_df[-1])
rownames(features_checked) = features_checked_df[,1]
# https://stackoverflow.com/questions/5738773/r-how-to-merge-two-matrix-according-to-their-column-and-row-names


# step 4a filter on outcome NA's

patientID_outcome = names(outcome_raw[complete.cases(outcome_raw)])
patientID_features = rownames(features_checked[complete.cases(features_checked),])
patientID_complete = intersect(patientID_outcome, patientID_features)

# patientID_outcome = names(outcome_raw)
# patientID_features = rownames(features_checked)
# patientID_complete = intersect(patientID_outcome, patientID_features)

# step 5a syntax input for prediction

outcome_complete = outcome_raw[patientID_complete]
patientID_complete = patientID_complete
features_complete = features_checked[patientID_complete,]
sum(clinical_raw$`Multiplex general panel` == 'yes', na.rm = TRUE) #89 patients multiplex
length(patientID_complete) # 60 patients viable for prediction

# outcome_complete = outcome_raw[patientID_complete]
# patientID_complete = patientID_complete
# features_complete = features_checked[patientID_complete,]



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
View(output)
saveRDS(object = output, file = 'prediction_tumors_borders_uni.RDS')




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
# [1] 0.4890286
# [1] 0.1682838
# 2.5% 
# 0.2 
# 97.5% 
# 0.8285714 




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

# groupings per image type
mygrouping1 = list(1:1366,1367:2732) 

out = grouping_prediction(mygrouping1)
aucs = out[[1]]
gammas = out[[2]]
mean(aucs); sd(aucs); quantile(aucs,0.025); quantile(aucs,0.975); #cross validated AUC
# [1] 0.6931429
# [1] 0.1524816
# 2.5% 
# 0.3335714 
# 97.5% 
# 0.9078571 
colMeans(gammas)/sum(colMeans(gammas)) #estimated group weights
# [1] 0.6507198 0.3492802

# groupings per groups of statistics general (distinct groups)
counts = which(str_detect(colnames(features),'counts_[a-z]+|density|X2stat'))  # L64 Regex: which(str_detect(colnames(features),'counts_[a-z]+|density|X2stat')) 
MED_MAD_SS = which(str_detect(colnames(features),'MED|MAD|distance_ratio')) # 108 Regex: which(str_detect(colnames(features),'MED|MAD|distance_ratio'))
centered = which(str_detect(colnames(features), 'Centered_')) # L1280 Regex: which(str_detect(colnames(features), 'Centered_'))
normalized = which(str_detect(colnames(features), 'Normalized_')) # L1280 Regex: which(str_detect(colnames(features), 'Normalized_'))

mygrouping2 = list(counts,MED_MAD_SS,normalized,centered) # L2732 sum(unlist(lapply(mygrouping1,function(x) length(x))))
# all features used: setdiff(colnames(features),unique(colnames(features)[unlist(mygrouping2)])) = empty

out = grouping_prediction(mygrouping2)
aucs = out[[1]]
gammas = out[[2]]
mean(aucs); sd(aucs); quantile(aucs,0.025); quantile(aucs,0.975); #cross validated AUC
# [1] 0.7445714
# [1] 0.1324375
# 2.5% 
# 0.4985714 
# 97.5% 
# 0.9078571 
colMeans(gammas)/sum(colMeans(gammas)) #estimated group weights
# [1] 0.37229102 0.54250764 0.06293505 0.02226629


# groupings per statistics detailed (distinct groups)
counts_ss = which(str_detect(colnames(features),'counts_[a-z]+|density|X2stat|distance_ratio'))  # L68 Regex: which(str_detect(colnames(features),'counts_[a-z]+|density|X2stat|distance_ratio')) 
MED = which(str_detect(colnames(features),'MED')) # L52 Regex: which(str_detect(colnames(features),'MED'))
# to view: View(colnames(features)[which(str_detect(colnames(features),'MED'))])
MAD = which(str_detect(colnames(features),'MAD')) # L52 Regex: which(str_detect(colnames(features),'MAD'))
Fstat = which(str_detect(colnames(features),'F_radius')) # L128 which(str_detect(colnames(features),'F_radius'))
Gstat = which(str_detect(colnames(features),'G_radius|Gdot_radius')) # L640 which(str_detect(colnames(features),'G_radius|Gdot_radius'))
Kstat = which(str_detect(colnames(features),'K_radius|Kdot_radius'))  # L640 which(str_detect(colnames(features),'K_radius|Kdot_radius'))
Lstat = which(str_detect(colnames(features),'L_radius|Ldot_radius')) # L640 which(str_detect(colnames(features),'L_radius|Ldot_radius'))
pcfstat = which(str_detect(colnames(features),'pcf_radius')) # L512 which(str_detect(colnames(features),'pcf_radius'))

mygrouping3 = list(counts_ss,MED,MAD,Fstat,Gstat,Kstat,Lstat,pcfstat) # L2732 sum(unlist(lapply(mygrouping,function(x) length(x))))
# all features used: setdiff(colnames(features),unique(colnames(features)[unlist(mygrouping3)])) = empty

out = grouping_prediction(mygrouping3)
aucs = out[[1]]
gammas = out[[2]]
mean(aucs); sd(aucs); quantile(aucs,0.025); quantile(aucs,0.975); #cross validated AUC
# [1] 0.744
# [1] 0.1404414
# 2.5% 
# 0.435 
# 97.5% 
# 0.9428571 
colMeans(gammas)/sum(colMeans(gammas)) #estimated group weights
# [1] 0.041109744 0.691394048 0.035005564 0.175659482 0.004272093 0.004496926 0.015256164 0.032805979


# groupings per statistics detailed with radius close (<10) and far (>=10) (distinct groups)
counts_ss = which(str_detect(colnames(features),'counts_[a-z]+|density|X2stat|distance_ratio')) # L68 Regex: which(str_detect(colnames(features),'counts_[a-z]+|density|X2stat|distance_ratio')) 
MED = which(str_detect(colnames(features),'MED')) # L52 Regex: which(str_detect(colnames(features),'MED'))
MAD = which(str_detect(colnames(features),'MAD')) # L52 Regex: which(str_detect(colnames(features),'MAD'))
Fstat_close = which(str_detect(colnames(features),'F_radius [0-9]_')) # L48
Fstat_far = which(str_detect(colnames(features),'F_radius [0-9][0-9]_')) # L80 
Gstat_close = which(str_detect(colnames(features),'G_radius [0-9]_|Gdot_radius [0-9]_')) # L240 
Gstat_far = which(str_detect(colnames(features),'G_radius [0-9][0-9]_|Gdot_radius [0-9][0-9]_')) # L400 
Kstat_close = which(str_detect(colnames(features),'K_radius [0-9]_|Kdot_radius [0-9]_')) # L240 
Kstat_far = which(str_detect(colnames(features),'K_radius [0-9][0-9]_|Kdot_radius [0-9][0-9]_')) # L400 
Lstat_close = which(str_detect(colnames(features),'L_radius [0-9]_|Ldot_radius [0-9]_')) # L240 
Lstat_far = which(str_detect(colnames(features),'L_radius [0-9][0-9]_|Ldot_radius [0-9][0-9]_')) # L400 
pcfstat_close = which(str_detect(colnames(features),'pcf_radius [0-9]_')) # L192 
pcfstat_far = which(str_detect(colnames(features),'pcf_radius [0-9][0-9]_')) # L320 

mygrouping4 = list(counts_ss,MED,MAD,
                  Fstat_close,Fstat_far,Gstat_close, Gstat_far,
                  Kstat_close,Kstat_far,Lstat_close,Lstat_far,
                  pcfstat_close,pcfstat_far) # L1366 sum(unlist(lapply(mygrouping,function(x) length(x))))
# all features used: setdiff(colnames(features),unique(colnames(features)[unlist(mygrouping)])) = empty

out = grouping_prediction(mygrouping4)
aucs = out[[1]]
gammas = out[[2]]
mean(aucs); sd(aucs); quantile(aucs,0.025); quantile(aucs,0.975); #cross validated AUC
# [1] 0.732
# [1] 0.1362177
# 2.5% 
# 0.435 
# 97.5% 
# 0.9428571 
colMeans(gammas)/sum(colMeans(gammas)) #estimated group weights
# [1] 0.0734123292 0.4080702416 0.0421463307 0.0048541466 0.2656108075 0.0000000000 0.0605658414 0.0260265917
# [9] 0.0136201945 0.0558714913 0.0005398924 0.0323981469 0.0168839862


# groupings per phenotype (overlapping groups)
Macrophage = which(str_detect(colnames(features),'Macrophage|distance_ratio')) #L1094
Tcells = which(str_detect(colnames(features),'Tcells|distance_ratio')) #L1094
Tumors = which(str_detect(colnames(features),'Tumor|distance_ratio')) #L1094
Others = which(str_detect(colnames(features),'Others')) #L1094

mygrouping5 = list(Macrophage,Tcells,Tumors,Others) # L5191 sum(unlist(lapply(mygrouping,function(x) length(x))))
# all features used: setdiff(colnames(features),unique(colnames(features)[unlist(mygrouping5)])) = empty

out = grouping_prediction(mygrouping5)
aucs = out[[1]]
gammas = out[[2]]
mean(aucs); sd(aucs); quantile(aucs,0.025); quantile(aucs,0.975); #cross validated AUC
# error
# error
# error
# error
# error
# error
colMeans(gammas)/sum(colMeans(gammas)) #estimated group weights
# error



## RANDOM FOREST WITH GROUP STRUCTURE (CORF)

