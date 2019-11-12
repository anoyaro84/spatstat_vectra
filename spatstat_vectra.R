

library(spatstat)
library(tidyverse)
library(phenoptr)
library(spdep)
library(zoo)

# phenotype can be converted by the "phenotype argument".
# why statistics was used as argument?
do_analyse <- function(Intable, PhenoOrder = NULL, Cols = NULL, phenotype = NULL, plotter = c(FALSE,FALSE), 
                       XposCol = 'Cell X Position', YposCol = 'Cell Y Position',
                       PhenoCol = 'Phenotype', sample_name = NULL, r_vec = NULL, ...) {
  
  csd <- Intable[, c(PhenoCol, XposCol, YposCol)]
  colnames(csd) = c('Phenotype', 'Cell X Position',  'Cell Y Position')
  pheno_vector = unique(csd$Phenotype)
  # print(pheno_vector)
  if (is.null(sample_name)) {
	  sample_name = 'Input sample'
  }
  
  if (is.null(PhenoOrder)) {
    PhenoOrder = pheno_vector # if no order is set, just take the order from data
  }
  
  if (is.element("",pheno_vector)){
    pheno_vector = pheno_vector[pheno_vector != ""]
  }
  pheno_vector = pheno_vector[order(match(pheno_vector,PhenoOrder))]
  
  # print(pheno_vector)
  
  # Replace phenotype label by definition
  if (!is.null(phenotype)) {
    for(pheno in names(phenotype)) {
      csd$Phenotype[[pheno]] = phenotype[[pheno]]
    }
  }
  
  if (is.null(r_vec)){
    dim_scale = min(max(csd[[XposCol]], max(csd[[YposCol]])))
    r_vec = dim_scale*c(0.1,0.9)
    # r_vec currently does not work well when given NULL due to several automated settings for the grid of r in the settings of each option
    # ASSUMPTION
    
  }
  
  print(paste("r_vec is", r_vec))
  
  # generate pairwise distance matrix for csd
  pairwise_distances = distance_matrix(Intable)
  
  # normal statistics: MAD and MED
  
  output = getMAD(Intable, pairwise_distances, pheno_vector)
  MED_min = output[1]
  MED = output[2]
  MAD_min = output[3]
  MAD = output[4]
  
  # normal statistics: Counts and Density
  Area_sample = max(csd[[XposCol]])*max(csd[[YposCol]])
  output = getDensity(csd, pheno_vector, Area_sample)

  counts_sample = output[1]
  density_sample = output[2]
  
  # normal statistics: Area and maxnorm
  dim_n = length(PhenoOrder)
  dim_m = 2*dim_n + 2
  
  statistics = as.data.frame(matrix(rep(0,dim_n*dim_m),dim_n, dim_m))
  colnames(statistics) = c(paste("Area K", PhenoOrder),"Area Kdot", paste("Maxnorm K", PhenoOrder), "Maxnorm Kdot")
  rownames(statistics) = PhenoOrder
  
  
  # Replace "+" with "T" and "-" with "F" in phenotype vector and in the csd for correct functioning of extracting inbuild statistics.
  pheno_vector = str_replace_all(pheno_vector,"[+]","T")
  pheno_vector = str_replace_all(pheno_vector,"[-]","F")
  
  
  csd$Phenotype = str_replace_all(csd$Phenotype,"[+]","T")
  csd$Phenotype = str_replace_all(csd$Phenotype,"[-]","F")
  
  n = length(pheno_vector)  
  #sample_name = str_replace(sample_name,".im3","")
  
  
  csd_ppp = ppp(x=csd[[XposCol]], y=csd[[YposCol]],
                window = owin(c(0, max(csd[[XposCol]])),
                              c(0, max(csd[[YposCol]]))),
                marks = as.factor(csd[[PhenoCol]]))
  if (plotter[1] == TRUE){
    plot1 = plot(csd_ppp, cols = Cols,
                 main = paste("Coordinates of cells and their phenotype in sample", sample_name), pch = 20)
  }
  

  # define all inbuild statistics
  options = c("G","F", "J","Gdot", "Jdot", "K", "L", "pcf", "Kdot", "Ldot")

  
  values_options = list()
  
  all_types_options_sample_name = list()
  
  for (option in options){
    # print(option)
    all_types = alltypes(csd_ppp,fun = paste(option), dataname = sample_name, envelope = FALSE)
    # print(all_types$fns)    
    all_types_options_sample_name[[option]] = all_types
    if (plotter[2] == TRUE){
      plot(all_types)
    }
    print(option)
    output = interpolate_r(all_types, r_vec, option)
    
    # r_close_list = output[1]
    # statistic_close_list = output[2]
    stop("round 1")
  }
  
  return(c(distances, output))
}


# function for retreiving the Median Absolute Deviation distance for each phenotype, and the Median Absolute Deviation for each phenotype to the closest other phenotype
getMAD <- function(data_with_distance, pairwise_distances, pheno_vector){
  
  MED = matrix(NA , nrow = length(pheno_vector), ncol = length(pheno_vector))
  colnames(MED) = pheno_vector
  rownames(MED) = pheno_vector
  
  MAD = MED
  MED_min = MED
  MAD_min = MED
  
  for (from in pheno_vector){
    filter_from = data_with_distance %>% filter(`Phenotype` == from & `Phenotype` != "")
    IDs_from = filter_from$`Cell ID`

    for (to in pheno_vector){
      distances_min = filter_from[[paste("Distance to",to)]]

      MED_min[paste(from), paste(to)] = median(distances_min)
      MAD_min[paste(from), paste(to)] = mad(distances_min)
      
      filter_to = data %>% filter(`Phenotype` == to & `Phenotype` != "")
      IDs_to = filter_to$`Cell ID`
      
      pairwise_to_from = pairwise_distances[IDs_from,IDs_to]
      
      MED[paste(from), paste(to)] = median(pairwise_to_from)
      MAD[paste(from), paste(to)] = mad(pairwise_to_from)
    }
  }
  return(c(MED_min, MED, MAD_min, MAD ))
}

# function for receiving the counts of the phenotypes in the sample and their density
getDensity <- function(data, pheno_vector, Area_sample){
  counts_sample = rep(0,length(pheno_vector))
  names(counts_sample) = pheno_vector
  
  density_sample = rep(0,length(pheno_vector))
  names(density_sample) = pheno_vector
  
  for (i in pheno_vector){
    n = dim(data %>% filter(`Phenotype` == i & `Phenotype` != ""))[1]
    counts_sample[[i]] = n
    density_sample[[i]] = n / Area_sample
  }
  
  return(c(counts_sample, density_sample))
}

# function interpolate statistic for given r
interpolate_r <- function(all_types, r_vec, option){
  r_close_list = list()
  statistic_close_list = list()

  for (r_i in r_vec){
    r_close_list[[option]][[paste("radius", r_i)]] = list()
    statistic_close_list[[option]][[paste("radius", r_i)]] = list()
    
    for (range in (1:length(all_types$fns))){
      r_close_list[[option]][[paste("radius", r_i)]][[paste(option, "fns which",range)]] = r_i
      
      
      if ((option == "K")|| (option == "L") || (option == "Kdot")|| (option == "Ldot")){
        
        stat1 = all_types[["fns"]][[range]][["border"]][left]
        stat2 = all_types[["fns"]][[range]][["border"]][right]
        print(c(stat1,stat2))
        
        a = (stat2-stat1)/(r2-r1) # y=ax+b
        b = stat2 - a*r2
        stat = a*(r_i-r2)+stat2
        statistic_close_list[[option]][[paste("radius", r_i)]][[option]][[paste(option, "fns which",range)]] = stat
      } else if (option == "pcf"){
        stat1 = all_types[["fns"]][[range]][["trans"]][left]
        stat2 = all_types[["fns"]][[range]][["trans"]][right]
        print(c(stat1,stat2))
        
        a = (stat2-stat1)/(r2-r1)
        b = stat2 - a*r2
        stat = a*(r_i-r2)+stat2
        
        statistic_close_list[[option]][[paste("radius", r_i)]][[option]][[paste(option, "fns which",range)]] = stat
      } else{
        
        view(all_types[["fns"]][[range]])
        r_emperic = all_types[["fns"]][[range]][["r"]]

        stat_emperic = all_types[["fns"]][[range]][["theohaz"]] #ASSUMPTION temporarily

        
        if (max(r_emperic)< r_i){
          warning("default r interval is to small")
          next(paste("skip", range, option))
        }
        
        dif = r_emperic - r_i
        dif_abs = abs(dif)
        condition = match(min(dif_abs),dif_abs)
        
        if (dif[condition]>0){
          left = condition - 1
          right = condition
        }else{
          left = condition
          right = condition + 1
        }
        
        r1 = r_emperic[left]
        r2 = r_emperic[right]
        
        print(c(r1,r2))
        stat1 = all_types[["fns"]][[range]][["theohaz"]][left]
        stat2 = all_types[["fns"]][[range]][["theohaz"]][right]
        
        # stat1 = all_types[["fns"]][[range]][["rs"]][left]
        # stat2 = all_types[["fns"]][[range]][["rs"]][right]
        print(c(stat1,stat2))
        
        a = (stat2-stat1)/(r2-r1)
        b = stat2 - a*r2
        stat = a*(r_i-r2)+stat2
        statistic_close_list[[option]][[paste("radius", r_i)]][[option]][[paste(option, "fns which",range)]] = stat
        
        print(paste("r left of r_i is", r1, ",r_i is", r_i, ",right of r_i is", r2))
        print(paste("stat left of r_i is", stat1, ",stat is", stat, ",right of stat is", stat2))
        stop("end")
      }
    }
  }
  return(c(r_close_list, statistic_close_list))
}


# Calculate the mean shortest distance for each pair of types, including itself
calculate_mean <- function(csd, pheno_vector, pheno_vector_absoluut, plotter){
  
  
  # calculate for each type of cell the mean distance to another type (inclusive itself). for all possible types in all the data.
  
  statistic_mean_sample_name = as.data.frame(matrix(rep(0,length(pheno_vector_absoluut)^2),c(length(pheno_vector_absoluut),length(pheno_vector_absoluut))))
  colnames(statistic_mean_sample_name) = paste("Mean shortest distance to",pheno_vector_absoluut)
  rownames(statistic_mean_sample_name) = pheno_vector_absoluut
  
  for (i in pheno_vector){
    if ("PAX5+PDL1+" %in% pheno_vector){
      dummy = csd %>%
        filter((`Phenotype` == i))
      statistic_mean_sample_name[i,"Mean shortest distance to PAX5+PDL1+"] = mean(dummy$`Distance to PAX5+PDL1+`)
      
      if (plotter[1] == TRUE){
        hist(csd$`Distance to PAX5+PDL1+`[csd$Phenotype == i],
             main = paste(i,"to PAX5+PDL1+"), ylab = "distance to PAX5+PDL1+")
        boxplot(csd$`Distance to PAX5+PDL1+`[csd$Phenotype == i],
                main = paste(i,"to PAX5+PDL1+"), ylab = "distance to PAX5+PDL1+")
      }
    }
    
    if ("PAX5+PDL1-" %in% pheno_vector){
      dummy = csd %>%
        filter((`Phenotype` == i))
      statistic_mean_sample_name[i,"Mean shortest distance to PAX5+PDL1-"] = mean(dummy$`Distance to PAX5+PDL1-`)
      
      if (plotter[1] == TRUE){
        hist(csd$`Distance to PAX5+PDL1-`[csd$Phenotype == i],
             main = paste(i,"to PAX5+PDL1-"), ylab = "distance to PAX5+PDL1-")
        boxplot(csd$`Distance to PAX5+PDL1-`[csd$Phenotype == i],
                main = paste(i,"to PAX5+PDL1-"), ylab = "distance to PAX5+PDL1-")
      }
    }
    
    if ("CD163+PDL1+" %in% pheno_vector){
      dummy = csd %>%
        filter((`Phenotype` == i))
      statistic_mean_sample_name[i,"Mean shortest distance to CD163+PDL1+"] = mean(dummy$`Distance to CD163+PDL1+`)
      
      if (plotter[1] == TRUE){
        hist(csd$`Distance to CD163+PDL1+`[csd$Phenotype == i],
             main = paste(i,"to CD163+PDL1+"), ylab = "distance to CD163+PDL1+")
        boxplot(csd$`Distance to CD163+PDL1+`[csd$Phenotype == i],
                main = paste(i,"to CD163+PDL1+"), ylab = "distance to CD163+PDL1+")
      }
    }
    
    if ("CD163+PDL1-" %in% pheno_vector){
      dummy = csd %>%
        filter((`Phenotype` == i))
      statistic_mean_sample_name[i,"Mean shortest distance to CD163+PDL1-"] = mean(dummy$`Distance to CD163+PDL1-`)
      
      if (plotter[1] == TRUE){
        hist(csd$`Distance to CD163+PDL1-`[csd$Phenotype == i],
             main = paste(i,"to CD163+PDL1-"), ylab = "distance to CD163+PDL1-")
        boxplot(csd$`Distance to CD163+PDL1-`[csd$Phenotype == i],
                main = paste(i,"to CD163+PDL1-"), ylab = "distance to CD163+PDL1-")
      }
    }
    
    if ("Other PDL1+" %in% pheno_vector){
      dummy = csd %>%
        filter((`Phenotype` == i))
      statistic_mean_sample_name[i,"Mean shortest distance to Other PDL1+"] = mean(dummy$`Distance to Other PDL1+`)
      
      if (plotter[1] == TRUE){
        hist(csd$`Distance to Other PDL1+`[csd$Phenotype == i],
             main = paste(i,"to Other PDL1+"), ylab = "distance to Other PDL1+")
        boxplot(csd$`Distance to Other PDL1+`[csd$Phenotype == i],
                main = paste(i,"to Other PDL1+"), ylab = "distance to Other PDL1+")
      }
    }
    
    if ("Other" %in% pheno_vector){
      dummy = csd %>%
        filter((`Phenotype` == i))
      statistic_mean_sample_name[i,"Mean shortest distance to Other"] = mean(dummy$`Distance to Other`)
      
      if (plotter[1] == TRUE){
        hist(csd$`Distance to Other`[csd$Phenotype == i],
             main = paste(i,"to Other"), ylab = "distance to Other")
        boxplot(csd$`Distance to Other`[csd$Phenotype == i],
                main = paste(i,"to Other"), ylab = "distance to Other")
      }
    }
    
    if ("CD3+CD8+PD1+" %in% pheno_vector){
      dummy = csd %>%
        filter((`Phenotype` == i))
      statistic_mean_sample_name[i,"Mean shortest distance to CD3+CD8+PD1+"] = mean(dummy$`Distance to CD163+PDL1+`)
      
      if (plotter[1] == TRUE){
        hist(csd$`Distance to CD3+CD8+PD1+`[csd$Phenotype == i])
        boxplot(csd$`Distance to CD3+CD8+PD1+`[csd$Phenotype == i],
                main = paste(i,"to  CD3+CD8+PD1+"), ylab = "distance to CD3+CD8+PD1+")
      }
    }
    
    if ("CD3+CD8+PD1-" %in% pheno_vector){
      dummy = csd %>%
        filter((`Phenotype` == i))
      statistic_mean_sample_name[i,"Mean shortest distance to CD3+CD8+PD1-"] = mean(dummy$`Distance to CD3+CD8+PD1-`)
      
      if (plotter[1] == TRUE){
        hist(csd$`Distance to CD3+CD8+PD1-`[csd$Phenotype == i],
             main = paste(i,"to  CD3+CD8+PD1-"), ylab = "distance to CD3+CD8+PD1-")
        boxplot(csd$`Distance to CD3+CD8+PD1-`[csd$Phenotype == i],
                main = paste(i,"to  CD3+CD8+PD1-"), ylab = "distance to CD3+CD8+PD1-")
      }
    }
    
    if ("CD3+CD8-PD1+" %in% pheno_vector){
      dummy = csd %>%
        filter((`Phenotype` == i))
      statistic_mean_sample_name[i,"Mean shortest distance to CD3+CD8-PD1+"] = mean(dummy$`Distance to CD3+CD8-PD1+`)
      
      if (plotter[1] == TRUE){
        hist(csd$`Distance to CD3+CD8-PD1+`[csd$Phenotype == i],
             main = paste(i,"to  CD3+CD8-PD1+"), ylab = "distance to CD3+CD8-PD1+")
        boxplot(csd$`Distance to CD3+CD8-PD1+`[csd$Phenotype == i],
                main = paste(i,"to  CD3+CD8-PD1+"), ylab = "distance to CD3+CD8-PD1+")
      }
    }
    
    if ("CD3+CD8-PD1-" %in% pheno_vector){
      dummy = csd %>%
        filter((`Phenotype` == i))
      statistic_mean_sample_name[i,"Mean shortest distance to CD3+CD8-PD1-"] = mean(dummy$`Distance to CD3+CD8-PD1-`)
      
      if (plotter[1] == TRUE){
        hist(csd$`Distance to CD3+CD8-PD1-`[csd$Phenotype == i],
             main = paste(i,"to  CD3+CD8-PD1-"), ylab = "distance to CD3+CD8-PD1-")
        boxplot(csd$`Distance to CD3+CD8-PD1-`[csd$Phenotype == i],
                main = paste(i,"to  CD3+CD8-PD1-"), ylab = "distance to CD3+CD8-PD1-")
      }
    }
  }
  return(statistic_mean_sample_name)
}

######### RUN function for calculating area and Maxnorm #####
calculate_area_norm <- function (ripleys, option){
  if (option == "K"){
    correction = "border"
  }
  if (option == "Kdot"){
    correction = "border"
  }
  x = ripleys[["r"]]
  y = abs(ripleys[["theo"]]-ripleys[[correction]])
  id = order(x)
  AUC = sum(diff(x[id])*rollmean(y[id],2))
  
  Norm_max = max(abs(ripleys[["theo"]]-ripleys[[correction]]))
  
  return(c(AUC,Norm_max))
}
