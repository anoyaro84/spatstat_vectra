#### load packages and github forks ####

library(tidyverse)
library(spatstat)
library(spdep)
library(remotes)
library(tiff)
library(phenoptr)
# library(zoo)
library(RColorBrewer)
library(reshape2)
library(latex2exp)



#### function master: do analyse on the path to the file ####

do_analyse <- function(seg_path, PhenoOrder = NULL, ColsOrder = NULL,
                       XposCol = 'Cell X Position', YposCol = 'Cell Y Position', PhenoCol = 'Phenotype',
                       sample_name = 'Input sample', plotter = c(FALSE,FALSE,FALSE), fig.prefix = '.',
                       r_vec = NULL, spatstat_statistics = NULL, ...) {
  
  
  # Create table with the right spatial dimensions such as described by the component file
  Intable = purrr::map_df(seg_path, read_cell_seg_data, pixels_per_micron = "auto",remove_units = FALSE)
  
  # replace empty phenotype with "Other"
  Intable$Phenotype[Intable$Phenotype == ""] = "Other"
  
  # define csd for ppp
  csd <- Intable[, c(PhenoCol, XposCol, YposCol)]
  colnames(csd) = c('Phenotype', 'Cell X Position',  'Cell Y Position')
  
  
  check_elsestate = FALSE
  
  if (is.null(PhenoOrder)) {
    PhenoOrder = unique(csd$Phenotype) # if no order is set, just take the order from data
    pheno_vector = PhenoOrder
    names(PhenoOrder) = PhenoOrder
    ColsOrder = brewer.pal(length(PhenoOrder), 'Set1') 
    
    colors_phenotype = ColsOrder
    names(colors_phenotype) = pheno_vector
  } else {
    if (!is.null(names(PhenoOrder))) {
      for (pheno in names(PhenoOrder)) {
        csd$Phenotype[csd$Phenotype %in% PhenoOrder[[pheno]]] = pheno
        Intable$Phenotype[Intable$Phenotype %in% PhenoOrder[[pheno]]] = pheno
      }
      pheno_vector = unique(csd$Phenotype)
    } else {
      PhenoOrder = unique(csd$Phenotype)
      pheno_vector = PhenoOrder
      names(PhenoOrder) = PhenoOrder
      ColsOrder = ColsOrder[names(ColsOrder) %in% PhenoOrder]
      
      check_elsestate = TRUE
    }
    colors_phenotype = ColsOrder
  }
  
  missing_in_data = setdiff(names(PhenoOrder),pheno_vector)
  if (!is_empty(missing_in_data)){
    warning('Target phenotype ', missing_in_data, ' is missing in the data sample\n')
  }
  
  if (TRUE %in% plotter) {
    output_dir <- file.path(fig.prefix, sample_name)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = T)
      cat("Directory created for", samplename, "with output directory", output_dir, fill = TRUE)
    } else {
      warning("Directory for ", samplename, " already exists with output directory ", output_dir, ". Figures were possibly overwritten.\n")
    }
  }
  
  if (is.null(r_vec)) {
    stop("Give at least one radius for parameter 'r_vec' to the function to start the analyse.\n")
  }
  
  
  spatstat_statistics_all = list("G","F", "J","Gdot", "Jdot", "K", "L", "pcf", "Kdot", "Ldot")
  
  if (is.null(spatstat_statistics)) {
    spatstat_statistics = spatstat_statistics_all
  } else if (all(spatstat_statistics %in% spatstat_statistics_all)) {
    spatstat_statistics = spatstat_statistics
  } else {
    stop("One or more spatial statistics in parameter 'spatstat_statistics' are not correctly defined")
  }
  
  
  
  #### normal statistics: Median and Median Absolute Deviation ####
  
  # Create Intable with nearest distances for each phenotype, here after the substitution of "" to "Other" and Simplyfying the Phenotypes
  Intable_with_distance = Intable %>%
      do(bind_cols(., find_nearest_distance(.)))
  cat("dimensions of the data with distances is ", dim(Intable_with_distance)[1], " times ", dim(Intable_with_distance)[2], fill = TRUE)
  
  # generate pairwise distance matrix for csd for use in getMAD
  pairwise_distance = distance_matrix(csd)
  
  # call getMAD function
  output = getMAD(Intable_with_distance, pairwise_distance, pheno_vector, missing_in_data)
  MED_min = output[[1]]
  MED = output[[2]]
  MAD_min = output[[3]]
  MAD = output[[4]]
  
  
  #### Creation Poisson Point Process and quadratcounts figures ####
  
  csd_ppp = ppp(x=csd[[XposCol]], y=csd[[YposCol]], 
                window = owin(c(min(csd[[XposCol]]), max(csd[[XposCol]])), c(min(csd[[YposCol]]), max(csd[[YposCol]]))),
                marks = factor(x = csd[[PhenoCol]], levels = pheno_vector)) #sort? pheno_vector[order(match(pheno_vector,names(PhenoOrder)))]
                # marks = factor(x = csd[[PhenoCol]], levels = names(PhenoOrder))) #sort? names(PhenoOrder)[order(match(names(PhenoOrder),pheno_vector))]
  unitname(csd_ppp) = list("micron", "microns", 1)
  
  
  if (isTRUE(plotter[[1]])) {
    
    png(filename = paste0(file.path(output_dir, sample_name),".png"), width = 600, height = 480)
    par(mar=rep(0.5, 4))
    plot(csd_ppp, cols = unlist(colors_phenotype[levels(csd_ppp$marks)]), xlab = "", ylab = "", main = "", pch = 20)
    title(paste("Location of cells and their phenotype\n in sample", sample_name), line = -3)
    dev.off()
  }
  
  
  ##### normal statistics: Counts and Density ####
  
  counts_sample = summary(csd_ppp)$marks[['frequency']]
  names(counts_sample) = levels(marks(csd_ppp))
  
  density_sample = summary(csd_ppp)$marks[['intensity']]
  names(density_sample) = levels(marks(csd_ppp))
  
  
  if (!is_empty(missing_in_data)){
    for (missing_pheno in missing_in_data){
      counts_sample[[missing_pheno]] = 0
      density_sample[[missing_pheno]] = 0
    }
  }
  # print(counts_sample)
  counts_normed_sample = counts_sample/csd_ppp$n
  names(counts_normed_sample) = names(counts_sample)
  
  
  phenos = c(pheno_vector, missing_in_data)
  # phenos = pheno_vector
  dim_square = length(phenos)
  
  counts_pairwise = matrix(NA , nrow = dim_square, ncol = dim_square)
  colnames(counts_pairwise) = phenos
  rownames(counts_pairwise) = phenos
  
  
  for (counter1 in seq_along(phenos)){
    phenotype1 = phenos[counter1]
    
    # use the symmetry of pairwise phenotypes and the inverse symmetry of the features for efficient loop
    sym_matrix_sequence = seq(counter1,seq_along(phenos)[dim_square])
    
    for (counter2 in sym_matrix_sequence){
      phenotype2 = phenos[counter2]
      
      if (phenotype1 %in% missing_in_data | phenotype2 %in% missing_in_data){
        counts_pairwise[phenotype1,phenotype2] = NA
      } else{
        if (phenotype1 == phenotype2){
          counts_pairwise[phenotype1,phenotype2] = 0 # symetric around 0
        } else {
          counts_pairwise[phenotype1,phenotype2] = log(counts_sample[[phenotype1]]/counts_sample[[phenotype2]])  # symetric around 0
          counts_pairwise[phenotype2,phenotype1] = log(counts_sample[[phenotype2]]/counts_sample[[phenotype1]])  # symetric around 0
        }
      }
    }
  }
  
  
  
  ##### normal statistics: Chi-squared statistics and quadratcount plot####
  
  quadratcount_X2statistic = list()
  quadratcount_X2statistic_normed = list()
  amount_pheno = length(pheno_vector)
  
  for (counter1 in seq_along(pheno_vector)){
    phenotype1 = pheno_vector[counter1]
    
    # use the symmetry of pairwise phenotypes for efficient loop
    sym_matrix_sequence = seq(counter1,seq_along(pheno_vector)[amount_pheno])
    
    for (counter2 in sym_matrix_sequence){
      phenotype2 = pheno_vector[counter2]
      
      # split ppp on pairwise phenotypes
      splitted = csd_ppp[(marks(csd_ppp) == phenotype1) | (marks(csd_ppp) == phenotype2)] 
      
      if (phenotype1 == phenotype2){
        # single phenotype Chi-squared quadratcount statistic and normalizing by the counts of corresponding phenotype
        quadrattest = quadrat.test(splitted) 
        quadratcount_X2statistic[[phenotype1]] = quadrattest$statistic[['X2']]
        quadratcount_X2statistic_normed[[phenotype1]] = quadrattest$statistic[['X2']]/counts_sample[[phenotype1]]
        
        if (isTRUE(plotter[[2]])){
          # plot quadratcounts for single phenotype and save in output directory
          png(filename = paste0(file.path(output_dir, sample_name),"_quadratcounts_", phenotype1, ".png"), width = 600, height = 480)
          par(mar=rep(0.5, 4))
          plot(splitted, cols = unlist(colors_phenotype[levels(csd_ppp$marks)]), xlab = "", ylab = "", main = "",  pch = 20)
          plot(quadratcount(splitted), add = TRUE)
          title(paste("Quadratcounts of", phenotype1, "\n in sample", sample_name), line = -3)
          dev.off()
          
          # plot single phenotype and save in output directory
          png(filename = paste0(file.path(output_dir, sample_name), '_', phenotype1, ".png"), width = 600, height = 480)
          par(mar=rep(0.5, 4))
          plot(splitted, cols = unlist(colors_phenotype[levels(csd_ppp$marks)]), xlab = "", ylab = "", main = "",  pch = 20)
          title(paste("Location of", phenotype1, "\n in sample", sample_name), line = -3)
          dev.off()
        }
      } else if (isTRUE(plotter[[2]])){
        # plot pairwise phenotypes and save in output directory
        png(filename = paste0(file.path(output_dir, sample_name), '_', phenotype1, '_', phenotype2, ".png"), width = 600, height = 480)
        par(mar=rep(0.5, 4)) # mar.panel=c(2,1,1,2)
        plot(splitted, cols = unlist(colors_phenotype[levels(csd_ppp$marks)]), xlab = "", ylab = "", main = "",  pch = 20)
        title(paste("Location of", phenotype1, "and", phenotype2, "\n in sample", sample_name), line = -3)
        dev.off()
      }
    }
  }
  
  if (!is_empty(missing_in_data)){
    for (missing_pheno in missing_in_data){
      quadratcount_X2statistic[[missing_pheno]] = NA
      quadratcount_X2statistic_normed[[missing_pheno]] = NA
    }
  }
  
  # browser()
  #### Replace "+" with "T" and "-" with "F" in all variables for correct functioning of extracting inbuild statistics ####
  
  PhenoOrder = lapply(PhenoOrder, function(x) {gsub("+", "T", x, fixed = TRUE)})
  PhenoOrder = lapply(PhenoOrder, function(x) {gsub("-", "F", x, fixed = TRUE)})
  if (isTRUE(check_elsestate)){
    names(PhenoOrder) = PhenoOrder
  }
  
  pheno_vector = gsub("+", "T", pheno_vector, fixed = TRUE)
  pheno_vector = gsub("-", "F", pheno_vector, fixed = TRUE)
  
  missing_in_data = gsub("+", "T", missing_in_data, fixed = TRUE)
  missing_in_data = gsub("-", "F", missing_in_data, fixed = TRUE)
  
  csd$Phenotype = sapply(csd$Phenotype, function(x) {gsub("+", "T", x, fixed = TRUE)})
  csd$Phenotype = sapply(csd$Phenotype, function(x) {gsub("-", "F", x, fixed = TRUE)})
  
  csd_ppp = ppp(x=csd[[XposCol]], y=csd[[YposCol]], 
                window = owin(c(min(csd[[XposCol]]), max(csd[[XposCol]])), c(min(csd[[YposCol]]), max(csd[[YposCol]]))),
                marks = factor(x = csd[[PhenoCol]], levels = pheno_vector)) #sort?
                # marks = factor(x = csd[[PhenoCol]], levels = names(PhenoOrder))) #sort?
  unitname(csd_ppp) = list("micron", "microns", 1)
  
  
  #### compute inbuild spatial statistics ####
  
  all_types_spatstat_statistics_sample_name = list()
  
  statistic_close_list = list()
  normalized_list = list()
  
  
  #### compute inbuild statistics with statistic-dependent correction-method ####
  for (spatstat_statistic in spatstat_statistics){
    
    cat('computing', spatstat_statistic, 'of alltypes', fill = TRUE)
    
    if (spatstat_statistic %in% list("K","L","Kdot","Ldot","pcf")){
      all_types = alltypes(csd_ppp,fun = paste(spatstat_statistic), correction = "iso", dataname = sample_name, envelope = TRUE, verb = TRUE)
    } else {
      all_types = alltypes(csd_ppp,fun = paste(spatstat_statistic), correction = "km", dataname = sample_name, envelope = TRUE, verb = TRUE)
    }
    
    # save object for debugging
    all_types_spatstat_statistics_sample_name[[spatstat_statistic]] = all_types
    
    # plot computation of inbuild statistic and save in output directory
    if (isTRUE(plotter[[3]])){
      
      png(filename = paste0(file.path(output_dir, sample_name),"_statistic_",spatstat_statistic,".png"), width = 720, height = 720)
      par(mar=rep(0.5, 4))
      plot(all_types, samex = TRUE)
      # plot(all_types)
      dev.off()
    }
    
    # interpolate the statistic value (and the normalized statistic value) for the user-defined radi in r_vec
    output = interpolate_r(all_types, r_vec, spatstat_statistic)
    
    statistic_close_list[[spatstat_statistic]] = output[[1]]
    normalized_list[[spatstat_statistic]] = output[[2]]
  }
  
  #### gather the output of the computations in a list ####
  output_data_raw = list()
  output_data_raw[["csd_ppp"]] = csd_ppp
  output_data_raw[["counts_sample"]] = counts_sample
  output_data_raw[["counts_normed_sample"]] = counts_normed_sample # new
  output_data_raw[["counts_pairwise"]] = counts_pairwise               # new
  # output_data_raw[["Area_sample"]] = Area_sample
  output_data_raw[["density_sample"]] = density_sample
  output_data_raw[["quadratcount_X2statistic"]] = quadratcount_X2statistic
  output_data_raw[["quadratcount_X2statistic_normed"]] = quadratcount_X2statistic_normed
  output_data_raw[["MED_min"]] = MED_min
  output_data_raw[["MED"]] = MED
  output_data_raw[["MAD_min"]] = MAD_min
  output_data_raw[["MAD"]] = MAD
  output_data_raw[["statistic_close_list"]] = statistic_close_list
  output_data_raw[["normalized_list"]] = normalized_list
  output_data_raw[["all_types_spatstat_statistics_sample_name"]] = all_types_spatstat_statistics_sample_name
  
  #### call feature_extract and output the prediction matrix ####
  
#  output_data_matrix = feature_extract(output_data_raw)
  
  #### gather both the raw data and the predection matrix data for the output data of do_analyse ####
  
#  output_data = list()
#  output_data[['output_data_raw']] = output_data_raw
#  output_data[['output_data_matrix']] = output_data_matrix
  
  return(output_data_raw)
}


#### function normal statistic: Median and Median Absolute Deviation ####
getMAD <- function(data_with_distance, pairwise_distances, pheno_vector, missing_in_data){
  
  phenos = c(pheno_vector, missing_in_data)
  # phenos = pheno_vector
  dim_square = length(phenos)
  
  MED = matrix(NA , nrow = dim_square, ncol = dim_square)
  colnames(MED) = phenos
  rownames(MED) = phenos
  
  MAD = MED
  MED_min = MED
  MAD_min = MED
  
  for (from in pheno_vector){
    filter_from = data_with_distance %>% filter(`Phenotype` == from)
    IDs_from = filter_from$`Cell ID`
    
    for (to in pheno_vector){
      distances_min = filter_from[[paste("Distance to",to)]]
      MED_min[paste(from), paste(to)] = median(distances_min)
      MAD_min[paste(from), paste(to)] = mad(distances_min)
      
      filter_to = data_with_distance %>% filter(`Phenotype` == to)
      IDs_to = filter_to$`Cell ID`
      
      pairwise_to_from = pairwise_distances[IDs_from,IDs_to]
      
      MED[paste(from), paste(to)] = median(pairwise_to_from)
      MAD[paste(from), paste(to)] = mad(pairwise_to_from)
    }
  }
  return(list(MED_min, MED, MAD_min, MAD))
}

# #### function normal statistic: Counts and density ####
# getDensity <- function(data, pheno_vector, Area_sample){
#   
#   counts_sample = rep(0,length(pheno_vector))
#   names(counts_sample) = pheno_vector
#   
#   density_sample = rep(0,length(pheno_vector))
#   names(density_sample) = pheno_vector
#   
#   for (phenotype in pheno_vector){
#     n = dim(data %>% filter(`Phenotype` == phenotype))[1]
#     counts_sample[[phenotype]] = n
#     density_sample[[phenotype]] = n / Area_sample
#   }
#   
#   return(list(counts_sample, density_sample))
# }

#### function interpolate spatial statistic: interpolate the (normalized) statistic value for the user-defined radi in r_vec ####
interpolate_r <- function(all_types, r_vec, spatstat_statistic){
  
  statistic_close_list = list()
  normalized_list = list()
  cat('interpolating',spatstat_statistic, fill = TRUE)
  
  # loop over every radius
  for (r_i in r_vec){
    
    statistic_close_list[[paste("radius", r_i)]] = list()
    normalized_list[[paste("radius", r_i)]] = list()
    
    # loop over every pairwise combination of phenotypes for the statistic
    for (index_pairwise in seq_along(all_types$fns)){
      
      statistic_pairwise_phenotypes =  all_types[["fns"]][[index_pairwise]]
      r_emperic = statistic_pairwise_phenotypes[["r"]]
      
      # browser()
      
      # find the surrounding r values for the user defined radius
      dif = r_emperic - r_i
      dif_abs = abs(dif)
      condition = match(min(dif_abs),dif_abs)
      
      if (dif[condition] > 0){
        left = condition - 1
        right = condition
      } else{
        left = condition
        right = condition + 1
      }
      # print(max(r_emperic))
      if (max(r_emperic)< r_i){
        # print('inside')
        cat('computed r interval (rmax = ', max(r_emperic),') is too small for user defined radius, inserted NA\n', fill = TRUE)
        
        statistic_close_list[[paste("radius", r_i)]][[paste(spatstat_statistic, "fns which",index_pairwise)]] = NA
        normalized_list[[paste("radius", r_i)]][[paste(spatstat_statistic, "fns which",index_pairwise)]] = NA
        
        next
      }
      # print('after check max')
      
      r1 = r_emperic[left]
      r2 = r_emperic[right]
      
      higher_bound = statistic_pairwise_phenotypes[["hi"]]
      high1 = higher_bound[left]
      high2 = higher_bound[right]
      
      a = (high2-high1)/(r2-r1) # y=ax+b
      b = high2 - a*r2
      high = a*(r_i-r2)+high2
      
      lower_bound = statistic_pairwise_phenotypes[["lo"]]
      low1 = lower_bound[left]
      low2 = lower_bound[right]
      
      a = (low2-low1)/(r2-r1) # y=ax+b
      b = low2 - a*r2
      low = a*(r_i-r2)+low2
      
      stat_theoretic = statistic_pairwise_phenotypes[["theo"]]
      stat_theo1 = stat_theoretic[left]
      stat_theo2 = stat_theoretic[right]
      a = (stat_theo2-stat_theo1)/(r2-r1) # y=ax+b
      b = stat_theo2 - a*r2
      stat_theo = a*(r_i-r2)+stat_theo2
      
      stat1 = statistic_pairwise_phenotypes[["obs"]][left]
      stat2 = statistic_pairwise_phenotypes[["obs"]][right]
      a = (stat2-stat1)/(r2-r1) # y=ax+b
      b = stat2 - a*r2
      stat = a*(r_i-r2)+stat2
      
      
      if(is.na(high) | is.na(low)){
        warning("NA in calculating significance bands, put in NA for normalized\n")
        normalized = NA
      }else if (high - low == 0){
        warning("dividing by 0 in normalizing, put in 0 for normalized\n")
        normalized = NA
      } else{
        normalized = (stat-stat_theo)/abs(high-low)
      }
      
      statistic_close_list[[paste("radius", r_i)]][[paste(spatstat_statistic, "fns which",index_pairwise)]] = stat - stat_theo
      normalized_list[[paste("radius", r_i)]][[paste(spatstat_statistic, "fns which",index_pairwise)]] = normalized
    }
  }
  # cat('done interpolating',spatstat_statistic, fill = TRUE)
  return(list(statistic_close_list, normalized_list))
}

#### function features: feature extract function ####
feature_extract <- function(outputs){
  
  cat('begin feature extraction', fill = TRUE)
  
  
  functions = c()
  rs = c()
  # get function names and rs
  for (out in outputs) {
    functions = union(functions, names(out$statistic_close_list))
    rs = union(rs, names(out$statistic_close_list[[1]]))
  }
  
  # get feature names
  feat_names = list()
  for (func in functions) {
    feat_names[[func]] = c()
    for (out in outputs) {
      feat_names[[func]] = union(feat_names[[func]],
                                 apply(expand.grid(dimnames(out$all_types_spatstat_statistics_sample_name[[func]]$which)), 
                                       1, function(x) gsub('/$', '', paste0(x, collapse='/')))
      )
    }
  }
  
  # create a matrix
  allfeat = lapply(feat_names, function(x) expand.grid(rs, x))
  allfeat = lapply(allfeat, function(x) {apply(x, 1, function(y) paste0(y, collapse='_'))})
  allfeat_flat = c()
  
  for (i in seq_along(allfeat)) {
    allfeat_flat = c(allfeat_flat, paste0(names(allfeat)[[i]], '_', allfeat[[i]]),
                     paste0('Normalized_',names(allfeat)[[i]], '_', allfeat[[i]])
    )
  }
  
  
  mat_ripleys = matrix(NA, nrow = length(allfeat_flat), ncol = length(outputs),
               dimnames = list(sort(allfeat_flat), names(outputs)))
 
  
  # fill matrix
  for (i in seq_along(outputs)) {
    out = outputs[[i]]
    name = names(outputs)[i]
    for (func in names(out$statistic_close_list)) {
      df = melt(out$all_types_spatstat_statistics_sample_name[[func]]$which)
      df$featname = gsub("/NA$", "", paste0(df$Var1, "/", df$Var2))
      for (r in rs) {
        data = as.data.frame(t(as.data.frame(out$statistic_close_list[[func]][[r]])))
        if (nrow(data) > 0){
          data$which = unlist(lapply(rownames(data), 
                                     function(x) as.numeric(tail(strsplit(x, '.', fixed=T)[[1]],1))))
          #ind = match(data$which, df$value)
          ind = match(df$value, data$which)
          df$Ffeatname = paste0(func, '_', r, '_', df$featname)
          df$measure = data$V1[ind]
          df$Nmeasure = as.data.frame(t(as.data.frame(out$normalized_list[[func]][[r]])))$V1[ind]
          mat_ripleys[df$Ffeatname, name] = df$measure
          mat_ripleys[paste0('Normalized_', df$Ffeatname), name] = df$Nmeasure
        }
      }
    }
  }
  
  
  # get phenotypes for counts and densities
  counts = c()
  counts_normed = c()
  dens = c()
  X2stat = c()
  X2stat_normed = c()
  
  for (out in outputs) {
    counts = union(counts, names(out$counts_sample))
    counts_normed = union(counts_normed,names(out$counts_normed_sample))
    dens = union(dens, names(out$density_sample))
    X2stat = union(X2stat, names(out$quadratcount_X2statistic))
    X2stat_normed = union(X2stat_normed, names(out$quadratcount_X2statistic_normed))
  }
  
  # create a matrix for counts
  mat_counts = matrix(0, nrow = length(counts), ncol = length(outputs),
               dimnames = list(paste0('counts_sample_', sort(counts)), names(outputs)))
  
  mat_counts_normed = matrix(0, nrow = length(counts_normed), ncol = length(outputs),
                      dimnames = list(paste0('counts_normed_sample_', sort(counts_normed)), names(outputs)))
  
  # create a matrix for density
  mat_density = matrix(0, nrow = length(dens), ncol = length(outputs),
                dimnames = list(paste0('density_sample_', sort(dens)), names(outputs)))
  
  # create a matrix for Chi-squared statistic of quadratcounts
  mat_X2stat = matrix(NA, nrow = length(X2stat), ncol = length(outputs),
                       dimnames = list(paste0('X2stat_sample_', sort(X2stat)), names(outputs)))
  
  # create a matrix for Chi-squared normlized statistic of quadratcounts
  mat_X2stat_normed = matrix(NA, nrow = length(X2stat_normed), ncol = length(outputs),
                      dimnames = list(paste0('X2stat_normed_sample_', sort(X2stat_normed)), names(outputs)))
  
  # fill matrices for counts and density
  for (i in seq_along(outputs)) {
    out = outputs[[i]]
    name = names(outputs)[i]
    data_counts = out$counts_sample
    data_counts_normed = out$counts_normed_sample
    data_density = out$density_sample
    data_X2stat = out$quadratcount_X2statistic
    data_X2stat_normed = out$quadratcount_X2statistic_normed
    for (featname in names(data_counts)){
      mat_counts[paste0('counts_sample_',featname),name] = data_counts[[featname]]
      mat_counts_normed[paste0('counts_normed_sample_',featname),name] = data_counts_normed[[featname]]
      mat_density[paste0('density_sample_',featname),name] = data_density[[featname]]
      mat_X2stat[paste0('X2stat_sample_',featname),name] = data_X2stat[[featname]]
      mat_X2stat_normed[paste0('X2stat_normed_sample_',featname),name] = data_X2stat_normed[[featname]]
    }
  }
  
  # pairwise counts
  
  phenos = c()
  
  for (out in outputs) {
    phenos = union(phenos, rownames(out$counts_pairwise))
  }
  
  counts_pairwise = lapply(phenos, function(x) expand.grid(x,phenos))
  counts_pairwise = lapply(counts_pairwise, function(x) {apply(x, 1, function(y) paste0(y, collapse='_'))})
  counts_pairwise_flat = c()
  
  for (i in seq_along(counts_pairwise)) {
    # counts_pairwise_flat = c(counts_pairwise_flat, paste0(names(counts_pairwise)[[i]], '_', counts_pairwise[[i]]))
    counts_pairwise_flat = c(counts_pairwise_flat, paste0(names(counts_pairwise)[[i]], counts_pairwise[[i]])) # why names(counts_pairwise)[[i]] ?
  }
  counts_pairwise = counts_pairwise_flat
  
  counts_pairwise_allfeat = paste0('counts_pairwise_', counts_pairwise)
  
  # create a matrix for the count of phenotype per phenotype for the features
  mat_counts_pairwise = matrix(NA, nrow = length(counts_pairwise_allfeat), ncol = length(outputs),
                   dimnames = list(counts_pairwise_allfeat, names(outputs)))
  
  # fill matrices
  for (i in seq_along(outputs)) {
    out = outputs[[i]]
    name = names(outputs)[i]
    data_counts_pairwise = out$counts_pairwise
    
    for (featname_from in rownames(data_counts_pairwise)){
      for (featname_to in colnames(data_counts_pairwise)){
        mat_counts_pairwise[paste0('counts_pairwise_', featname_from, '_', featname_to),name] = data_counts_pairwise[featname_from, featname_to]
      }
    }
  }
  
  
  
  
  # get phenotypes for median minimal, median, MAD minimal, MAD.
  MED_min_pheno = c()
  MED_pheno = c()
  MAD_min_pheno = c()
  MAD_pheno = c()
  
  for (out in outputs) {
    MED_min_pheno = union(MED_min_pheno, rownames(out$MED_min))
    MED_pheno = union(MED_pheno, rownames(out$MED))
    MAD_min_pheno = union(MAD_min_pheno, rownames(out$MAD_min))
    MAD_pheno = union(MAD_pheno, rownames(out$MAD))
  }
  
  collect_pheno = sort(union(MED_min_pheno, union(MED_pheno, union(MAD_min_pheno,MAD_pheno))))
  
  
  allfeat_min = lapply(collect_pheno, function(x) expand.grid(collect_pheno, x))
  allfeat_min = lapply(allfeat_min, function(x) {apply(x, 1, function(y) paste0(y, collapse='_'))})
  allfeat_min_flat = c()
  
  for (i in seq_along(allfeat_min)) {
    # allfeat_min_flat = c(allfeat_min_flat, paste0(names(allfeat_min)[[i]], '_', allfeat_min[[i]]))
    allfeat_min_flat = c(allfeat_min_flat, paste0(names(allfeat_min)[[i]], allfeat_min[[i]])) # why names(allfeat_min)[[i]] ?
  }
  allfeat_min = allfeat_min_flat
  
  allfeat_normal = combn(collect_pheno,2, simplify = FALSE)
  allfeat_normal_identity = lapply(collect_pheno,function(x) rep(x,2))
  allfeat_normal = c(allfeat_normal, allfeat_normal_identity)
  allfeat_normal = sort(sapply(allfeat_normal, function(y) paste0(y, collapse = '_')))
  
  
  MED_min_allfeat = paste0('MED_min_', allfeat_min)
  MED_allfeat = paste0('MED_', allfeat_normal)
  MAD_min_allfeat = paste0('MAD_min_', allfeat_min)
  MAD_allfeat = paste0('MAD_', allfeat_normal)
  
  
  # create a matrix for MED_min NON-symmetric for the features
  mat_med_min = matrix(NA, nrow = length(MED_min_allfeat), ncol = length(outputs),
                dimnames = list(MED_min_allfeat, names(outputs)))
  
  # create a matrix for MED symmetric for the features
  mat_med = matrix(NA, nrow = length(MED_allfeat), ncol = length(outputs),
                dimnames = list(MED_allfeat, names(outputs)))
  
  # create a matrix for MAD_min NON-symmetric for the features
  mat_mad_min = matrix(NA, nrow = length(MAD_min_allfeat), ncol = length(outputs),
                dimnames = list(MAD_min_allfeat, names(outputs)))
  
  # create a matrix for MAD symmetric for the features
  mat_mad = matrix(NA, nrow = length(MAD_allfeat), ncol = length(outputs),
                dimnames = list(MAD_allfeat, names(outputs)))
  
  
  # fill matrices
  for (i in seq_along(outputs)) {
    out = outputs[[i]]
    name = names(outputs)[i]
    data_MED_min = out$MED_min
    
    data_MED = out$MED
    data_MED = data_MED[sort(colnames(data_MED)),sort(rownames(data_MED))]
    data_MAD_min = out$MAD_min
    data_MAD = out$MAD
    data_MAD = data_MAD[sort(colnames(data_MAD)),sort(rownames(data_MAD))]
    
    for (featname_from in rownames(data_MED_min)){
      for (featname_to in colnames(data_MED_min)){
        mat_med_min[paste0('MED_min_', featname_from, '_', featname_to),name] = data_MED_min[featname_from, featname_to]
      }
    }
    
    for (row in seq_along(rownames(data_MED))){
      featname_from = rownames(data_MED)[row]
      for (col in row:length(colnames(data_MED))){
        featname_to = colnames(data_MED)[col]
        mat_med[paste0('MED_', featname_from, '_', featname_to),name] = data_MED[featname_from, featname_to]
      }
    }
    
    for (featname_from in rownames(data_MAD_min)){
      for (featname_to in colnames(data_MAD_min)){
        mat_mad_min[paste0('MAD_min_', featname_from, '_', featname_to),name] = data_MAD_min[featname_from, featname_to]
      }
    }
    
    for (row in seq_along(rownames(data_MAD))){
      featname_from = rownames(data_MAD)[row]
      for (col in row:length(colnames(data_MAD))){
        featname_to = colnames(data_MAD)[col]
        mat_mad[paste0('MAD_', featname_from, '_', featname_to),name] = data_MAD[featname_from, featname_to]
      }
    }
  }
  
  mat = t(rbind(mat_counts,mat_counts_normed, mat_counts_pairwise, mat_density, mat_X2stat, mat_X2stat_normed,
              mat_med_min, mat_med, mat_mad_min, mat_mad, mat_ripleys))
  
  cat('end feature extraction')
  
  return(mat)
}


#### NOT RUN Calculate the mean shortest distance for each pair of types, including itself ####
# calculate_mean <- function(csd, pheno_vector, pheno_vector_absoluut, plotter){
#     
#     
#     # calculate for each type of cell the mean distance to another type (inclusive itself). for all possible types in all the data.
#     
#     statistic_mean_sample_name = as.data.frame(matrix(rep(0,length(pheno_vector_absoluut)^2),c(length(pheno_vector_absoluut),length(pheno_vector_absoluut))))
#     colnames(statistic_mean_sample_name) = paste("Mean shortest distance to",pheno_vector_absoluut)
#     rownames(statistic_mean_sample_name) = pheno_vector_absoluut
# 
#     for (i in pheno_vector){
#         if ("PAX5+PDL1+" %in% pheno_vector){
#             dummy = csd %>%
#             filter((`Phenotype` == i))
#             statistic_mean_sample_name[i,"Mean shortest distance to PAX5+PDL1+"] = mean(dummy$`Distance to PAX5+PDL1+`)
# 
#             if (plotter[1] == TRUE){
#                 hist(csd$`Distance to PAX5+PDL1+`[csd$Phenotype == i],main = paste(i,"to PAX5+PDL1+"), ylab = "distance to PAX5+PDL1+")
#                 boxplot(csd$`Distance to PAX5+PDL1+`[csd$Phenotype == i],main = paste(i,"to PAX5+PDL1+"), ylab = "distance to PAX5+PDL1+")
#             }
#         }
#         
#         if ("PAX5+PDL1-" %in% pheno_vector){
#             dummy = csd %>%
#             filter((`Phenotype` == i))
#             statistic_mean_sample_name[i,"Mean shortest distance to PAX5+PDL1-"] = mean(dummy$`Distance to PAX5+PDL1-`)
# 
#             if (plotter[1] == TRUE){
#                 hist(csd$`Distance to PAX5+PDL1-`[csd$Phenotype == i],main = paste(i,"to PAX5+PDL1-"), ylab = "distance to PAX5+PDL1-")
#                 boxplot(csd$`Distance to PAX5+PDL1-`[csd$Phenotype == i],main = paste(i,"to PAX5+PDL1-"), ylab = "distance to PAX5+PDL1-")
#             }
#         }
#         
#         if ("CD163+PDL1+" %in% pheno_vector){
#           dummy = csd %>%
#             filter((`Phenotype` == i))
#           statistic_mean_sample_name[i,"Mean shortest distance to CD163+PDL1+"] = mean(dummy$`Distance to CD163+PDL1+`)
#           
#           if (plotter[1] == TRUE){
#             hist(csd$`Distance to CD163+PDL1+`[csd$Phenotype == i],
#                  main = paste(i,"to CD163+PDL1+"), ylab = "distance to CD163+PDL1+")
#             boxplot(csd$`Distance to CD163+PDL1+`[csd$Phenotype == i],
#                     main = paste(i,"to CD163+PDL1+"), ylab = "distance to CD163+PDL1+")
#           }
#         }
#         
#         if ("CD163+PDL1-" %in% pheno_vector){
#           dummy = csd %>%
#             filter((`Phenotype` == i))
#           statistic_mean_sample_name[i,"Mean shortest distance to CD163+PDL1-"] = mean(dummy$`Distance to CD163+PDL1-`)
#           
#           if (plotter[1] == TRUE){
#             hist(csd$`Distance to CD164+PDL1-`[csd$Phenotype == i],
#                  main = paste(i,"to CD163+PDL1-"), ylab = "distance to CD163+PDL1-")
#             boxplot(csd$`Distance to CD163+PDL1-`[csd$Phenotype == i],
#                     main = paste(i,"to CD163+PDL1-"), ylab = "distance to CD163+PDL1-")
#           }
#         }
#         
#         if ("Other PDL1+" %in% pheno_vector){
#           dummy = csd %>%
#             filter((`Phenotype` == i))
#           statistic_mean_sample_name[i,"Mean shortest distance to Other PDL1+"] = mean(dummy$`Distance to Other PDL1+`)
#           
#           if (plotter[1] == TRUE){
#             hist(csd$`Distance to Other PDL1+`[csd$Phenotype == i],
#                  main = paste(i,"to Other PDL1+"), ylab = "distance to Other PDL1+")
#             boxplot(csd$`Distance to Other PDL1+`[csd$Phenotype == i],
#                     main = paste(i,"to Other PDL1+"), ylab = "distance to Other PDL1+")
#           }
#         }
#         
#         if ("Other" %in% pheno_vector){
#           dummy = csd %>%
#             filter((`Phenotype` == i))
#           statistic_mean_sample_name[i,"Mean shortest distance to Other"] = mean(dummy$`Distance to Other`)
#           
#           if (plotter[1] == TRUE){
#             hist(csd$`Distance to Other`[csd$Phenotype == i],
#                  main = paste(i,"to Other"), ylab = "distance to Other")
#             boxplot(csd$`Distance to Other`[csd$Phenotype == i],
#                     main = paste(i,"to Other"), ylab = "distance to Other")
#           }
#         }
#         
#         if ("CD3+CD8+PD1+" %in% pheno_vector){
#           dummy = csd %>%
#             filter((`Phenotype` == i))
#           statistic_mean_sample_name[i,"Mean shortest distance to CD3+CD8+PD1+"] = mean(dummy$`Distance to CD163+PDL1+`)
#           
#           if (plotter[1] == TRUE){
#             hist(csd$`Distance to CD3+CD8+PD1+`[csd$Phenotype == i])
#             boxplot(csd$`Distance to CD3+CD8+PD1+`[csd$Phenotype == i],
#                     main = paste(i,"to  CD3+CD8+PD1+"), ylab = "distance to CD3+CD8+PD1+")
#           }
#         }
#         
#         if ("CD3+CD8+PD1-" %in% pheno_vector){
#           dummy = csd %>%
#             filter((`Phenotype` == i))
#           statistic_mean_sample_name[i,"Mean shortest distance to CD3+CD8+PD1-"] = mean(dummy$`Distance to CD3+CD8+PD1-`)
#           
#           if (plotter[1] == TRUE){
#             hist(csd$`Distance to CD3+CD8+PD1-`[csd$Phenotype == i],
#                  main = paste(i,"to  CD3+CD8+PD1-"), ylab = "distance to CD3+CD8+PD1-")
#             boxplot(csd$`Distance to CD3+CD8+PD1-`[csd$Phenotype == i],
#                     main = paste(i,"to  CD3+CD8+PD1-"), ylab = "distance to CD3+CD8+PD1-")
#           }
#         }
#         
#         if ("CD3+CD8-PD1+" %in% pheno_vector){
#           dummy = csd %>%
#             filter((`Phenotype` == i))
#           statistic_mean_sample_name[i,"Mean shortest distance to CD3+CD8-PD1+"] = mean(dummy$`Distance to CD3+CD8-PD1+`)
#           
#           if (plotter[1] == TRUE){
#             hist(csd$`Distance to CD3+CD8-PD1+`[csd$Phenotype == i],
#                  main = paste(i,"to  CD3+CD8-PD1+"), ylab = "distance to CD3+CD8-PD1+")
#             boxplot(csd$`Distance to CD3+CD8-PD1+`[csd$Phenotype == i],
#                     main = paste(i,"to  CD3+CD8-PD1+"), ylab = "distance to CD3+CD8-PD1+")
#           }
#         }
#         
#         if ("CD3+CD8-PD1-" %in% pheno_vector){
#           dummy = csd %>%
#             filter((`Phenotype` == i))
#           statistic_mean_sample_name[i,"Mean shortest distance to CD3+CD8-PD1-"] = mean(dummy$`Distance to CD3+CD8-PD1-`)
#           
#           if (plotter[1] == TRUE){
#             hist(csd$`Distance to CD3+CD8-PD1-`[csd$Phenotype == i],
#                  main = paste(i,"to  CD3+CD8-PD1-"), ylab = "distance to CD3+CD8-PD1-")
#             boxplot(csd$`Distance to CD3+CD8-PD1-`[csd$Phenotype == i],
#                     main = paste(i,"to  CD3+CD8-PD1-"), ylab = "distance to CD3+CD8-PD1-")
#           }
#         }
#     }
#     return(statistic_mean_sample_name)
# }

#### NOT RUN function for calculating area and Maxnorm ####
calculate_area_norm <- function (ripleys, spatstat_statistic){
  if (spatstat_statistic == "K"){
    correction = "border"
  }
  if (spatstat_statistic == "Kdot"){
    correction = "border"
  }
  x = ripleys[["r"]]
  y = abs(ripleys[["theo"]]-ripleys[[correction]])
  id = order(x)
  AUC = sum(diff(x[id])*rollmean(y[id],2))
  
  Norm_max = max(abs(ripleys[["theo"]]-ripleys[[correction]]))
  
  return(c(AUC,Norm_max))
}
