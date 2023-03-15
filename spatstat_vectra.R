#### load packages and github forks ####

library(tidyverse)
library(spatstat)
library(spdep)
#library(remotes)
#library(tiff)
library(phenoptr)
# library(zoo)
library(RColorBrewer)
library(reshape2)
library(latex2exp)



#### function master: do analyse on the path to the file ####

do_analyse <- function(seg_path, PhenoOrder = NULL, ColsOrder = NULL,
                       XposCol = 'Cell X Position', YposCol = 'Cell Y Position', PhenoCol = 'Phenotype',
                       sample_name = 'Input sample', plotter = c(FALSE,FALSE,FALSE), fig.prefix = '.',
                       r_vec = NULL, spatstat_statistics = 'ALL', 
                       reference = 'Tumors',plotOnly=NULL,image_shape = 'concave', ...) {
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
    print(colors_phenotype)
    print(PhenoOrder)
    names(colors_phenotype) =  PhenoOrder
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
  
  if (isTRUE(spatstat_statistics == 'ALL')) {
    spatstat_statistics = spatstat_statistics_all
  } else if (all(spatstat_statistics %in% spatstat_statistics_all)) {
    spatstat_statistics = spatstat_statistics
  } else if (is.null(spatstat_statistics)) {
    spatstat_statistics = list()
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
  output = getMAD(Intable_with_distance, pairwise_distance, pheno_vector, missing_in_data, reference=reference)
  MED_min = output[[1]]
  MED = output[[2]]
  MAD_min = output[[3]]
  MAD = output[[4]]
  ratio_distances = output[[5]]
  
  
  #### Creation Poisson Point Process and quadratcounts figures ####
    
    if(image_shape == 'circle'){
        window <- disc(radius = min(2*sd(csd[[XposCol]]), 2*sd(csd[[YposCol]])), 
                              centre = c(mean(range((csd[[XposCol]]))), mean(range(csd[[YposCol]]))))
    }else{
        # default shape: rectangle
        window <- owin(c(min(csd[[XposCol]]), max(csd[[XposCol]])), c(min(csd[[YposCol]]), max(csd[[YposCol]])))
    }

   if(!image_shape %in% c('convex','concave')){
  csd_ppp = ppp(x=csd[[XposCol]], y=csd[[YposCol]], 
                window = window,
                marks = factor(x = csd[[PhenoCol]], levels = pheno_vector)) #sort? pheno_vector[order(match(pheno_vector,names(PhenoOrder)))]
                                        # marks = factor(x = csd[[PhenoCol]], levels = names(PhenoOrder))) #sort? names(PhenoOrder)[order(match(names(PhenoOrder),pheno_vector))]
   }else if(image_shape == 'convex'){
       csd_ppp = ppp(x=csd[[XposCol]], y=csd[[YposCol]], 
                     poly = list(x=convexhull.xy(x=csd[[XposCol]], y=csd[[YposCol]])$bdry[[1]]$x,
                                 y=convexhull.xy(x=csd[[XposCol]], y=csd[[YposCol]])$bdry[[1]]$y),
                     marks = factor(x = csd[[PhenoCol]], levels = pheno_vector)) 
       

   }else if(image_shape == 'concave'){
       concave <- concaveman::concaveman(as.matrix(csd[,c(XposCol,YposCol)]),concavity = 0.685,length_threshold=50)
       csd_ppp = ppp(x=csd[[XposCol]], y=csd[[YposCol]],
                     poly = list(apply(concave,2, rev)),
                     marks = factor(x = csd[[PhenoCol]], levels = pheno_vector))
   }
    unitname(csd_ppp) = list("micron", "microns", 1)
    
    
    if (isTRUE(plotter[[1]])) {
        
        png(filename = paste0(file.path(output_dir, sample_name),".png"), width = 600, height = 480)
        par(mar=rep(0.5, 4))
        plot(csd_ppp, cols = unlist(colors_phenotype[levels(csd_ppp$marks)]), xlab = "", ylab = "", main = "", pch = 20)
        title(paste("Location of cells and their phenotype\n in sample", sample_name), line = -3)
        dev.off()
    }
    if(plotOnly != TRUE){
  ##### normal statistics: Counts and Density ####
  
  counts_sample = summary(csd_ppp)$marks[['frequency']]
  names(counts_sample) = levels(marks(csd_ppp))
  
  density_sample = summary(csd_ppp)$marks[['intensity']]
  names(density_sample) = levels(marks(csd_ppp))
  print(counts_sample)
  
  if (!is_empty(missing_in_data)){
    for (missing_pheno in missing_in_data){
      counts_sample[[missing_pheno]] = 0
      density_sample[[missing_pheno]] = 0
    }
  }
  
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
    for (counter2 in seq_along(phenos)){
      phenotype2 = phenos[counter2]
      
      if (phenotype2 %in% missing_in_data){
        counts_pairwise[phenotype1,phenotype2] = NA
      } else if (phenotype1 %in% missing_in_data){
        counts_pairwise[phenotype1,phenotype2] = 0
      } else{
        counts_pairwise[phenotype1,phenotype2] = counts_sample[[phenotype1]]/counts_sample[[phenotype2]]
        #counts_pairwise[phenotype1,phenotype2] = log2(counts_sample[[phenotype1]]+1)-log2(counts_sample[[phenotype2]]+1)
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
  if(!is.null(spatstat_statistics)){
    for (spatstat_statistic in spatstat_statistics){
      
      cat('computing', spatstat_statistic, 'of alltypes', fill = TRUE)
      
      if (spatstat_statistic %in% list("K","L","Kdot","Ldot","pcf")){
        all_types = alltypes(csd_ppp,fun = paste(spatstat_statistic), envelope = TRUE, correction = "iso", dataname = sample_name, verb = FALSE,reuse = FALSE)
      } else {
        all_types = alltypes(csd_ppp,fun = paste(spatstat_statistic), envelope = TRUE, correction = "km", dataname = sample_name, verb = FALSE,reuse = FALSE)
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
  }
  #### gather the output of the computations in a list ####
  output_data_raw = list()
  output_data_raw[["csd"]] = csd
  output_data_raw[["distance_matrix"]] = pairwise_distance 
  output_data_raw[["csd_ppp"]] = csd_ppp
  output_data_raw[["counts_sample"]] = counts_sample
  output_data_raw[["counts_normed_sample"]] = counts_normed_sample
  output_data_raw[["counts_pairwise"]] = counts_pairwise
  output_data_raw[["density_sample"]] = density_sample
  output_data_raw[["quadratcount_X2statistic"]] = quadratcount_X2statistic
  output_data_raw[["quadratcount_X2statistic_normed"]] = quadratcount_X2statistic_normed
  output_data_raw[["MED_min"]] = MED_min
  output_data_raw[["MED"]] = MED
  output_data_raw[["MAD_min"]] = MAD_min
  output_data_raw[["MAD"]] = MAD
  output_data_raw[["ratio_distances"]] = ratio_distances
  
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
}

#### function normal statistic: Median and Median Absolute Deviation ####
getMAD <- function(data_with_distance, pairwise_distances, pheno_vector, missing_in_data,
                   reference = "Tumors" # used for relative distance calculation
                   ){
  
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
      
      MED[paste(from), paste(to)] = median(pairwise_to_from[pairwise_to_from > 0]) # median(pairwise_to_from[pairwise_to_from > 0]
      MAD[paste(from), paste(to)] = mad(pairwise_to_from[pairwise_to_from > 0])
    }
  }
  

  ratio_distances = NULL 
  combination = combn(setdiff(pheno_vector, reference), 2)
  combination = cbind(combination,   rbind(combination[2,], combination[1,]))

  

  if (!any(reference %in% missing_in_data)) {
      refs = data_with_distance %>% filter(`Phenotype` %in% reference)

      if (nrow(refs)!=0) {
          for (i in 1:ncol(combination)) {
              if ((!any(combination[,i] %in% missing_in_data))) {
                  Ctype1 = combination[1,i]
                  Ctype2 = combination[2,i]
                  dist_ctype1 = data_with_distance %>% filter(`Phenotype` == Ctype1)

         
                  mat = matrix(NA, nrow=length(refs$`Cell ID`), ncol=5,
                           dimnames = list(paste0('Cell ID ', refs$`Cell ID`),
                                           c(paste0('Cell ID ',Ctype1), paste0('Distance to ', Ctype1),
                                             paste0('Cell ID ',Ctype2), paste0('Distance to ', Ctype2),
                                             paste0('Relative_distance_', Ctype1, '_', Ctype2)
                                           )
                                    )
                               )
                  for (refID in refs$`Cell ID`) {
                      ctype1_ID = as.integer(refs[which(refs$`Cell ID` == refID), paste0('Cell ID ', Ctype1)])
                      mindist_ctype1 = as.numeric(refs[which(refs$`Cell ID` == refID), paste0('Distance to ', Ctype1)])                  
                      ctype2_ID = as.integer(dist_ctype1[which(dist_ctype1$`Cell ID` == ctype1_ID), 
                                         paste0('Cell ID ', Ctype2)])

                      mindist_ctype1_ctype2 = as.numeric(dist_ctype1[which(dist_ctype1$`Cell ID` == ctype1_ID), 
                                                  paste0('Distance to ', Ctype2)])
                      mat[paste0('Cell ID ', refID), ] =
                          c(ctype1_ID, mindist_ctype1, ctype2_ID, mindist_ctype1_ctype2,
                            mindist_ctype1/mindist_ctype1_ctype2
                            )
                  }
                  ratio_distances = cbind(ratio_distances, mat)
              }
          }
      }
  }


  #ratio_distances = matrix(NA, nrow = 1, ncol = 5, 
  #                         dimnames = list(paste('Cell ID'), c('Distance Tumor to Tcell','Cell ID Tcell','Distance Tcell to Macrophage', 'Cell ID Macrophage', 'distance_ratio_Tumor_Tcell_Macrophage')))
  
  #if (!any('Tumors' %in% missing_in_data)){
  #  tumors = data_with_distance %>% filter(`Phenotype` == 'Tumors')
    
  #  ratio_distances = matrix(NA, nrow = length(tumors$`Cell ID`), ncol = 5, 
  #                           dimnames = list(paste('Cell ID', tumors$`Cell ID`),c('Distance Tumor to Tcell','Cell ID Tcell','Distance Tcell to Macrophage', 'Cell ID Macrophage', 'distance_ratio_Tumor_Tcell_Macrophage')))
    
  #  if ((!any(c('Tcells', 'Macrophage') %in% missing_in_data))){
  #    tcells = data_with_distance %>% filter(`Phenotype` == 'Tcells')
  #    macrophages = data_with_distance %>% filter(`Phenotype` == 'Macrophage')
  #    
  #    for (tumor_ID in tumors$`Cell ID`){
  #      tcell_ID = as.integer(tumors[which(tumors$`Cell ID` == tumor_ID), 'Cell ID Tcells'])
  #      ratio_distances[paste('Cell ID', tumor_ID),'Cell ID Tcell'] = tcell_ID
  #      distance_tumor_to_tcell = as.numeric(tumors[which(tumors$`Cell ID` == tumor_ID), 'Distance to Tcells'])
  #      ratio_distances[paste('Cell ID', tumor_ID),'Distance Tumor to Tcell'] = distance_tumor_to_tcell
        
       # macrophage_ID = as.integer(tcells[which(tcells$`Cell ID` == tcell_ID), 'Cell ID Macrophage'])
  #      ratio_distances[paste('Cell ID', tumor_ID),'Cell ID Macrophage'] = macrophage_ID
  #      distance_tcell_to_macrophage = as.numeric(tcells[which(tcells$`Cell ID` == tcell_ID), 'Distance to Macrophage'])
  #      ratio_distances[paste('Cell ID', tumor_ID),'Distance Tcell to Macrophage'] = distance_tcell_to_macrophage
  #      
  #      ratio_distances[paste('Cell ID', tumor_ID),'distance_ratio_Tumor_Tcell_Macrophage'] = distance_tumor_to_tcell/distance_tcell_to_macrophage # spatial Score SS
  #    }
  #  }
  #}
  
  
  return(list(MED_min, MED, MAD_min, MAD, ratio_distances))
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
      statistic_pairwise_phenotypes =  as.fv(all_types[["fns"]][[index_pairwise]]) 
      # Different hardcoding maybe in future, more intuitive on which pair is selected: double loop (row,col) over phenotypes in all_types.
      # statistic_pairwise_phenotypes =  as.fv(all_types[row,col])
      # has implications for definition statistic_close_list, normalized_list and feature_extract-function. 
      
      # View(statistic_pairwise_phenotypes)
      # browser()
      ftheo = as.function(statistic_pairwise_phenotypes, value = 'theo', extrapolate = TRUE)
      stat_theo = ftheo(r_i)
      
      if (any(is.na(statistic_pairwise_phenotypes$obs[-1]))){
        # NA or NaN in 'obs': too few counts in observed pattern to compute centered and normalized statistic
        cat('NA or NaN \'obs\': too few counts in observed pattern to compute centered and normalized statistic, both set to NA',fill = T)
        centered = NA
        normalized = NA
      } else {
        # enough counts in observed pattern to compute centered statistic
        fobs = as.function(statistic_pairwise_phenotypes, value = 'obs', extrapolate = TRUE)
        stat_obs = fobs(r_i)
        centered = stat_obs-stat_theo
        if (any(is.na(statistic_pairwise_phenotypes$lo[-1]) | is.na(statistic_pairwise_phenotypes$hi[-1]))){
          # NA or NaN in 'lo' and/or 'hi': too few counts in observed pattern to compute normalized statistic
          cat("NA in calculating significance bands so width significance band is infinit, normalized statistic set to 0.", fill = T)
          normalized = 0
        } else {
          # enough counts in observed pattern to compute significance band, and thus centered statistic
          fenv = as.function(statistic_pairwise_phenotypes, value = c('lo','hi'), extrapolate = TRUE)
          low = fenv(r_i,'lo')
          high = fenv(r_i,'hi')
          if (isTRUE(high == low)){
            width_eps = 10^(-5)
            cat('significance band width is 0, instead normalizing by width_eps = 10^(-5).', fill = T)
            normalized = centered/width_eps
          } else {
            width = abs(high-low)
            normalized = centered/width
          }
        }
      }
      
      statistic_close_list[[paste("radius", r_i)]][[paste(spatstat_statistic, "fns which",index_pairwise)]] = centered
      normalized_list[[paste("radius", r_i)]][[paste(spatstat_statistic, "fns which",index_pairwise)]] = normalized
    }
  }
  return(list(statistic_close_list, normalized_list))
}

#### function features: feature extract function ####
feature_extract <- function(outputs){
  
  cat('begin feature extraction', fill = TRUE)
  #browser()
  spatstat_statistics_available = outputs[[1]][['statistic_close_list']]
  if (!is.null(spatstat_statistics_available)){
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
 
  mat_counts_lognormed = matrix(0, nrow = length(counts_normed), ncol = length(outputs),
                      dimnames = list(paste0('counts_lognormed_sample_', sort(counts_normed)), names(outputs)))
  
  # create a matrix for density
  mat_density = matrix(0, nrow = length(dens), ncol = length(outputs),
                dimnames = list(paste0('density_sample_', sort(dens)), names(outputs)))
 
  mat_density_lognormed = matrix(0, nrow = length(dens), ncol = length(outputs),
                dimnames = list(paste0('density_lognormed_sample_', sort(dens)), names(outputs)))
  
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
    size_image <- diff(out$csd_ppp$window$xrange) * diff(out$csd_ppp$window$yrange)
    data_counts = out$counts_sample/size_image
    data_counts_normed = out$counts_normed_sample
    data_counts_lognormed = log(out$counts_normed_sample*10^6+1)
    data_density = out$density_sample
    data_density_normed = log(out$density_sample*10^6+1)
    data_X2stat = out$quadratcount_X2statistic
    data_X2stat_normed = out$quadratcount_X2statistic_normed
    for (featname in names(data_counts)){
      mat_counts[paste0('counts_sample_',featname),name] = data_counts[[featname]]
      mat_counts_normed[paste0('counts_normed_sample_',featname),name] = data_counts_normed[[featname]]
      mat_counts_lognormed[paste0('counts_lognormed_sample_',featname),name] = data_counts_lognormed[[featname]]
      mat_density[paste0('density_sample_',featname),name] = data_density[[featname]]
      mat_density_lognormed[paste0('density_lognormed_sample_',featname),name] = data_density_normed[[featname]]
      mat_X2stat[paste0('X2stat_sample_',featname),name] = data_X2stat[[featname]]
      mat_X2stat_normed[paste0('X2stat_normed_sample_',featname),name] = data_X2stat_normed[[featname]]
    }
  }
  
  # get phenotypes for pairwise counts
  
  phenos = c()
  
  for (out in outputs) {
    phenos = union(phenos, rownames(out$counts_pairwise))
  }
  
  counts_pairwise = outer(X = phenos, Y = phenos, FUN = 'paste', sep = '_')
  
  counts_pairwise_flat = c()
  for (i in seq_along(phenos)){
    for (j in seq_along(phenos)){
      if (i != j){
        counts_pairwise_flat = c(counts_pairwise_flat, paste0(counts_pairwise[i,j]))
      }
    }
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
    data_counts_pairwise = log(out$counts_pairwise+1)
    
    for (featname_from in rownames(data_counts_pairwise)){
      for (featname_to in colnames(data_counts_pairwise)){
        if (featname_from != featname_to){
          mat_counts_pairwise[paste0('counts_pairwise_', featname_from, '_', featname_to),name] = data_counts_pairwise[featname_from, featname_to]
        }
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
  
  allfeat_min = outer(X = collect_pheno, Y = collect_pheno, FUN = 'paste', sep = '_')
  
  allfeat_min_flat = c()
  for (i in seq_along(collect_pheno)){
    for (j in seq_along(collect_pheno)){
      allfeat_min_flat = c(allfeat_min_flat, paste0(allfeat_min[i,j]))
    }
  }
  
  # allfeat_min = lapply(collect_pheno, function(x) expand.grid(collect_pheno, x))
  # allfeat_min = lapply(allfeat_min, function(x) {apply(x, 1, function(y) paste0(y, collapse='_'))})
  # allfeat_min_flat = c()
  # 
  # for (i in seq_along(allfeat_min)) {
  #   # allfeat_min_flat = c(allfeat_min_flat, paste0(names(allfeat_min)[[i]], '_', allfeat_min[[i]]))
  #   allfeat_min_flat = c(allfeat_min_flat, paste0(names(allfeat_min)[[i]], allfeat_min[[i]])) # why names(allfeat_min)[[i]] ?
  # }
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

  # Identify all possible pairs of cell types (for pairwise relative distances)
  CellPairs = c()
  for (out in outputs) {
      CellPairs = union(CellPairs, grep('Relative_distance_', colnames(out$ratio_distances), value=T))
  }

  # create a matrix for pairwise distances 
  mat_ratio_distances_median = matrix(NA, nrow = length(CellPairs), ncol = length(outputs),
                       dimnames = list(CellPairs, names(outputs)))
  
  # create a matrix for MAD symmetric for the features
  mat_ratio_distances_mad = matrix(NA, nrow = length(CellPairs), ncol = length(outputs),
                                   dimnames = list(CellPairs, names(outputs)))
  
  
  # fill matrices
  for (i in seq_along(outputs)) {
    out = outputs[[i]]
    name = names(outputs)[i]
    data_mat_ratio_distances = out$ratio_distances
    for (feat in CellPairs) {
        if (feat %in% colnames(data_mat_ratio_distances)){
            mat_ratio_distances = data_mat_ratio_distances[,feat]
            mat_ratio_distances_median[feat, name] = median(mat_ratio_distances)
            mat_ratio_distances_mad[feat, name] = mad(mat_ratio_distances)
        }
    }
  }
  
  
  if (!is.null(spatstat_statistics_available)){
    mat = t(rbind(mat_counts,mat_counts_normed, mat_counts_lognormed, mat_counts_pairwise, mat_density, mat_density_lognormed,  
                  mat_X2stat, mat_X2stat_normed,
                mat_med_min, mat_med, mat_mad_min, mat_mad, mat_ratio_distances_median, mat_ratio_distances_mad, mat_ripleys))
  } else {
    mat = t(rbind(mat_counts,mat_counts_normed,mat_counts_lognormed, mat_counts_pairwise, mat_density, mat_density_lognormed, 
                  mat_X2stat, mat_X2stat_normed,
                  mat_med_min, mat_med, mat_mad_min, mat_mad, mat_ratio_distances_median, mat_ratio_distances_mad))
  }
  
  cat('end feature extraction')
  
  return(mat)
}

#### RUN function statisticPerPatient: take the statistic (mean or median) over all the features for each patient ####

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

