

library(spatstat)
library(tidyverse)
library(phenoptr)
library(spdep)
library(zoo)

# phenotype can be converted by the "phenotype argument".
# why statistics was used as argument?
do_analyse <- function(Intable, PhenoOrder = NULL, Cols = NULL, phenotype = NULL, plotter = c(FALSE,FALSE,FALSE), 
                       fig.prefix = './', XposCol = 'Cell X Position', YposCol = 'Cell Y Position',
                       PhenoCol = 'Phenotype', sample_name = NULL, r_vec = NULL, envelope_bool = TRUE, ...) {
  
  csd <- Intable[, c(PhenoCol, XposCol, YposCol)]
  colnames(csd) = c('Phenotype', 'Cell X Position',  'Cell Y Position')
  pheno_vector = unique(csd$Phenotype)
  print(pheno_vector)
  
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
  
  
  # generate pairwise distance matrix for csd and filter csd on ""
  pairwise_distance_all = distance_matrix(csd)

  rows = select_rows(csd,pheno_vector)
  pairwise_distance_filtered = pairwise_distance_all[rows,]
  csd <- csd[rows, ]

  
  # normal statistics: MAD and MED
  # uses Intable as it needs the minimal distances in the table and their ID's
  
  output = getMAD(Intable, pairwise_distance_all, pheno_vector)
  MED_min = output[[1]]
  MED = output[[2]]
  MAD_min = output[[3]]
  MAD = output[[4]]
  
  # normal statistics: Counts and Density
  Area_sample = max(csd[[XposCol]])*max(csd[[YposCol]])
  output = getDensity(csd, pheno_vector, Area_sample)

  counts_sample = output[[1]]
  density_sample = output[[2]]
  
  
  # normal statistics: Area and maxnorm
  dim_n = length(PhenoOrder)
  dim_m = 2*dim_n + 2
  
  # statistics = as.data.frame(matrix(rep(0,dim_n*dim_m),dim_n, dim_m))
  # colnames(statistics) = c(paste("Area K", PhenoOrder),"Area Kdot", paste("Maxnorm K", PhenoOrder), "Maxnorm Kdot")
  # rownames(statistics) = PhenoOrder
  
  
  # Replace "+" with "T" and "-" with "F" in phenotype vector and in the csd for correct functioning of extracting inbuild statistics.
  pheno_vector = str_replace_all(pheno_vector,"[+]","T")
  pheno_vector = str_replace_all(pheno_vector,"[-]","F")
  
  
  csd$Phenotype = str_replace_all(csd$Phenotype,"[+]","T")
  csd$Phenotype = str_replace_all(csd$Phenotype,"[-]","F")
  
  n = length(pheno_vector)  
  
  
  
  csd_ppp = ppp(x=csd[[XposCol]], y=csd[[YposCol]],
                window = owin(c(0, max(csd[[XposCol]])),
                              c(0, max(csd[[YposCol]]))),
                marks = as.factor(csd[[PhenoCol]]))
  if (plotter[1] == TRUE){
    png(filename = paste0(fig.prefix, sample_name,".png"))
    par(mar = c(0,2,0,0)+0.1)
    plot(csd_ppp, cols = Cols, xlab = "", ylab = "", main = "", pch = 20)
    title(paste("Location of cells and their phenotype\n in sample", sample_name), line = -5)
    dev.off()
  }
  
  
  quadratcount_pvalue = list()

  for (phenotype in pheno_vector){
    splitted = csd_ppp[marks(csd_ppp) == phenotype]
    quadrattest = quadrat.test(splitted)
    quadratcount_pvalue[[phenotype]] = quadrattest$p.value
    
    if (plotter[2] == TRUE){
      png(filename = paste0(fig.prefix, sample_name,"_quadratcounts_",phenotype,".png"))
      par(mar = c(0,2,0,0)+0.1)
      # par(mar=c(6,6,6,6)+0.1, mgp = c(2,1,0))
      plot(splitted, cols = Cols, xlab = "", ylab = "", main = "",  pch = 20)
      plot(quadratcount(splitted), add = TRUE)
      title(paste("Counts of", phenotype, "in sample", sample_name), line = -5)
      dev.off()
      }
  }
  
  
  
  # define all inbuild statistics
  options = c("G","F", "J","Gdot", "Jdot", "K", "L", "pcf", "Kdot", "Ldot")
  # G and F are alreay in J.
  
  values_options = list()
  
  all_types_options_sample_name = list()
 
  r_close_list = list()
  statistic_close_list = list()
  normalized_list = list()
  
  
  for (option in options){
    
    if (option %in% c("K","L","Kdot","Ldot","pcf")){  # (option == "K")|| (option == "L") || (option == "Kdot")|| (option == "Ldot") || (option == "pcf")){
      all_types = alltypes(csd_ppp,fun = paste(option), dataname = sample_name, envelope = envelope_bool, correction = "iso")
    } else{
      all_types = alltypes(csd_ppp,fun = paste(option), dataname = sample_name, envelope = envelope_bool, correction = "km")
    }

    
    all_types_options_sample_name[[option]] = all_types
    
    if (plotter[3] == TRUE){
      if (option %in% c("G","J","K", "L","pcf")){ # (option == "K")|| (option == "L") || (option == "pcf")){
        width = 800
        height = 700
        res = 80
      } else{
        width = 800
        height = 900
        res = 80
      }
        # par(mar = c(1,2,1,1)+0.1, oma = c())
      # par(mar=c(10,10,10,10)+0.1, mgp = c(3,1,0))
      # par(mar = c(5,4,4,2) + 0.1, oma = c(1,1,1,1))
      png(filename = paste0(fig.prefix, sample_name,"_statistic_",option,".png"),
          width = width, height = height, res = res)
      plot(all_types)
      dev.off()
    }
    
    output = interpolate_r(all_types, r_vec, option, envelope_bool)
    
    r_close_list[[option]] = output[[1]]
    statistic_close_list[[option]] = output[[2]]
    normalized_list[[option]] = output[[3]]
    
  }
  
  output_all = list()
  output_all[["pairwise_distance_filtered"]] = pairwise_distance_filtered
  output_all[["counts_sample"]] = counts_sample
  output_all[["Area_sample"]] = Area_sample
  output_all[["density_sample"]] = density_sample
  output_all[["quadratcount_pvalue"]] = quadratcount_pvalue
  output_all[["MED_min"]] = MED_min
  output_all[["MED"]] = MED
  output_all[["MAD_min"]] = MAD_min
  output_all[["MAD"]] = MAD
  output_all[["r_close_list"]] = r_close_list
  output_all[["statistic_close_list"]] = statistic_close_list
  output_all[["normalized_list"]] = normalized_list
  output_all[["all_types_options_sample_name"]] = all_types_options_sample_name

  return(output_all)
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
  return(list(MED_min, MED, MAD_min, MAD))
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
  
  return(list(counts_sample, density_sample))
}

# function interpolate statistic for given r
interpolate_r <- function(all_types, r_vec, option, envelope_bool){
  
  r_close_list = list()
  statistic_close_list = list()
  normalized_list = list()
  
  print(option)
  
  for (r_i in r_vec){
    
    r_close_list[[paste("radius", r_i)]] = list()
    statistic_close_list[[paste("radius", r_i)]] = list()
    normalized_list[[paste("radius", r_i)]] = list()
    
    for (range in (1:length(all_types$fns))){
      
      
      r_close_list[[paste("radius", r_i)]][[paste(option, "fns which",range)]] = r_i
      
      r_emperic = all_types[["fns"]][[range]][["r"]]
      dif = r_emperic - r_i
      dif_abs = abs(dif)
      condition = match(min(dif_abs),dif_abs)
  
      if (dif[condition] > 0){
        left = condition - 1
        right = condition
      }else{
        left = condition
        right = condition + 1
      }

      r1 = r_emperic[left]
      r2 = r_emperic[right]

      if ((option == "K")|| (option == "L") || (option == "Kdot")|| (option == "Ldot") || (option == "pcf")){
        
        
        
        if (max(r_emperic)< r_i){
          warning("default r interval is to small")
          next(paste("skip", range, option))
        }
        
        
        if (envelope_bool == TRUE){
          
          
          
          higher_bound = all_types[["fns"]][[range]][["hi"]]
          high1 = higher_bound[left]
          high2 = higher_bound[right]
          # print(paste("high1",c(high1,high2)))
          a = (high2-high1)/(r2-r1) # y=ax+b
          b = high2 - a*r2
          high = a*(r_i-r2)+high2
          
          lower_bound = all_types[["fns"]][[range]][["lo"]]
          low1 = lower_bound[left]
          low2 = lower_bound[right]
          # print(paste("low1",c(low1,low2)))
          a = (low2-low1)/(r2-r1) # y=ax+b
          b = low2 - a*r2
          low = a*(r_i-r2)+low2
          
          stat_theoretic = all_types[["fns"]][[range]][["theo"]]
          stat_theo1 = stat_theoretic[left]
          stat_theo2 = stat_theoretic[right]
          a = (stat_theo2-stat_theo1)/(r2-r1) # y=ax+b
          b = stat_theo2 - a*r2
          stat_theo = a*(r_i-r2)+stat_theo2
          
          stat1 = all_types[["fns"]][[range]][["obs"]][left]
          stat2 = all_types[["fns"]][[range]][["obs"]][right]
          a = (stat2-stat1)/(r2-r1) # y=ax+b
          b = stat2 - a*r2
          stat = a*(r_i-r2)+stat2
          
        if(is.na(high) || is.na(low)){
          warning("NA in all_types[[\"fns\"]][[range]], put in NaN")
          normalized = NaN
        }else if (high - low == 0){
            warning("deviding by 0")
            normalized = NaN
        } else{
            normalized = (stat-stat_theo)/abs(high-low)
        }
          
          statistic_close_list[[paste("radius", r_i)]][[paste(option, "fns which",range)]] = stat
          normalized_list[[paste("radius", r_i)]][[paste(option, "fns which",range)]] = normalized
          
        } # else{ # envelope_bool == FALSE Not needed?
# 
#           browser()
# 
#           a = (stat_theo2-stat_theo1)/(r2-r1) # y=ax+b
#           b = stat_theo2 - a*r2
#           stat_theo = a*(r_i-r2)+stat_theo2
# 
#           stat1 = all_types[["fns"]][[range]][["iso"]][left]
#           stat2 = all_types[["fns"]][[range]][["iso"]][right]
# 
#           a = (stat2-stat1)/(r2-r1) # y=ax+b
#           b = stat2 - a*r2
#           stat = a*(r_i-r2)+stat2
#         }

        
      } else{ # option == "G","F", "J","Gdot", "Jdot"
        
        # browser()
        
        if (max(r_emperic)< r_i){
          warning("default r interval is to small")
          next(paste("skip", range, option))
        }
        
        if (envelope_bool == TRUE){
          
          
          # browser()
          higher_bound = all_types[["fns"]][[range]][["hi"]]
          high1 = higher_bound[left]
          high2 = higher_bound[right]
          # print(paste("high1",c(high1,high2)))
          a = (high2-high1)/(r2-r1) # y=ax+b
          b = high2 - a*r2
          high = a*(r_i-r2)+high2
          
          lower_bound = all_types[["fns"]][[range]][["lo"]]
          low1 = lower_bound[left]
          low2 = lower_bound[right]
          # print(paste("low1",c(low1,low2)))
          a = (low2-low1)/(r2-r1) # y=ax+b
          b = low2 - a*r2
          low = a*(r_i-r2)+low2
          
          stat_theoretic = all_types[["fns"]][[range]][["theo"]]
          stat_theo1 = stat_theoretic[left]
          stat_theo2 = stat_theoretic[right]
          a = (stat_theo2-stat_theo1)/(r2-r1) # y=ax+b
          b = stat_theo2 - a*r2
          stat_theo = a*(r_i-r2)+stat_theo2
          
          stat1 = all_types[["fns"]][[range]][["obs"]][left]
          stat2 = all_types[["fns"]][[range]][["obs"]][right]
          a = (stat2-stat1)/(r2-r1) # y=ax+b
          b = stat2 - a*r2
          stat = a*(r_i-r2)+stat2
          
          if(is.na(high) || is.na(low)){
            warning("NA in all_types[[\"fns\"]][[range]], put in NaN")
            normalized = NaN
          }else if (high - low == 0){
            warning("deviding by 0")
            normalized = NaN
          } else{
            normalized = (stat-stat_theo)/abs(high-low)
          }
          
          statistic_close_list[[paste("radius", r_i)]][[paste(option, "fns which",range)]] = stat
          normalized_list[[paste("radius", r_i)]][[paste(option, "fns which",range)]] = normalized
          
        } # else{ # envelope_bool == FALSE
          # 
          # browser()
          # 
          # stat1 = all_types[["fns"]][[range]][["km"]][left]
          # stat2 = all_types[["fns"]][[range]][["km"]][right]
          # 
          # a = (stat2-stat1)/(r2-r1) # y=ax+b
          # b = stat2 - a*r2
          # stat = a*(r_i-r2)+stat2
        # }
        # browser()
        # stat1 = all_types[["fns"]][[range]][["iso"]][left]
        # stat2 = all_types[["fns"]][[range]][["iso"]][right]
        # 
        # a = (stat2-stat1)/(r2-r1) # y=ax+b
        # b = stat2 - a*r2
        # stat = a*(r_i-r2)+stat2
        # 
        # statistic_close_list[[paste("radius", r_i)]][[paste(option, "fns which",range)]] = stat
        # 
        # # print(paste("r left of r_i is", r1, ",r_i is", r_i, ",right of r_i is", r2))
        # # print(paste("stat left of r_i is", stat1, ",stat is", stat, ",right of stat is", stat2))
        
        
      }
        # browser()
        # stat_emperic = all_types[["fns"]][[range]][["km"]] #ASSUMPTION temporarily
        # 
        # stat1 = all_types[["fns"]][[range]][["km"]][left]
        # stat2 = all_types[["fns"]][[range]][["km"]][right]
        # 
        # a = (stat2-stat1)/(r2-r1)
        # b = stat2 - a*r2
        # stat = a*(r_i-r2)+stat2
        # 
        # statistic_close_list[[paste("radius", r_i)]][[paste(option, "fns which",range)]] = stat
        # 
        # (stat-theo)/|high-low|
        
        
        # print(paste("r left of r_i is", r1, ",r_i is", r_i, ",right of r_i is", r2))
        # print(paste("stat left of r_i is", stat1, ",stat is", stat, ",right of stat is", stat2))
        
    }
  }
  
  return(list(r_close_list, statistic_close_list, normalized_list))
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
        hist(csd$`Distance to CD164+PDL1-`[csd$Phenotype == i],
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


feature_extract <- function(outputs){
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
					apply(expand.grid(dimnames(out$all_types_options_sample_name[[func]]$which)), 
					      1, function(x) gsub('/$', '', paste0(x, collapse='/')))
				)
		}
	}

	# create a matrix
	allfeat = lapply(feat_names, function(x) expand.grid(rs, x))
	allfeat = lapply(allfeat, function(x) {apply(x, 1, function(y) paste0(y, collapse='_'))})
	allfeat_flat = c()
	for (i in 1:length(allfeat)) {
		allfeat_flat = c(allfeat_flat, paste0(names(allfeat)[[i]], '_', allfeat[[i]]),
				 paste0('Normalized_',names(allfeat)[[i]], '_', allfeat[[i]])
				 )
	}


	mat = matrix(NA, nrow = length(allfeat_flat), ncol = length(outputs),
		     dimnames = list(sort(allfeat_flat), names(outputs)))

	# fill matrix
	for (i in 1:length(outputs)) {
		out = outputs[[i]]
		name = names(outputs)[i]
		for (func in names(out$statistic_close_list)) {
			df = melt(out$all_types_options_sample_name[[func]]$which)
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
					mat[df$Ffeatname, name] = df$measure
					mat[paste0('Normalized_', df$Ffeatname), name] = df$Nmeasure
				}
			}
		}
	}

	# add more features (density, count, etc..)
#	extra_features <- c('MED_min

	return(mat)
}


