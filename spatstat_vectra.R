

library(spatstat)
library(tidyverse)
library(phenoptr)
library(spdep)
library(zoo)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(latex2exp)

# phenotype can be converted by the "phenotype argument".
# why statistics was used as argument?
do_analyse <- function(Intable, PhenoOrder = NULL, ColsOrder = NULL,
                       XposCol = 'Cell X Position', YposCol = 'Cell Y Position', PhenoCol = 'Phenotype',
                       sample_name = NULL, plotter = c(FALSE,FALSE,FALSE), fig.prefix = '.', fig.width, fig.height,
                       r_vec = NULL, options = NULL, envelope_bool = TRUE, ...) {
    
    
    
    csd <- Intable[, c(PhenoCol, XposCol, YposCol)]
    colnames(csd) = c('Phenotype', 'Cell X Position',  'Cell Y Position')
    
    # replace empty phenotype with "Other"
    csd$Phenotype[csd$Phenotype == ""] = "Other"
    Intable$Phenotype[Intable$Phenotype == ""] = "Other"
    
    
    if (is.null(sample_name)) {
      sample_name = 'Input sample'
    }
    
    
    check_elsestate = FALSE
    
    if (is.null(PhenoOrder)) {
        PhenoOrder = unique(csd$Phenotype) # if no order is set, just take the order from data
        pheno_vector = PhenoOrder
        names(PhenoOrder) = PhenoOrder
        ColsOrder = brewer.pal(length(PhenoOrder), 'Set1') 
        
        colors_phenotype = ColsOrder
        names(colors_phenotype) = pheno_vector
    }else {
        if (!is.null(names(PhenoOrder))) {
            # for (pheno_replace in names(PhenoOrder)) {
            # 
            #     pheno = PhenoOrder[names(PhenoOrder) == pheno_replace]
            #     csd$Phenotype[csd$Phenotype %in% pheno]  = pheno_replace
            #     Intable$Phenotype[Intable$Phenotype %in% pheno] = pheno_replace
            # }
            for (pheno in names(PhenoOrder)) {
                csd$Phenotype[csd$Phenotype %in% PhenoOrder[[pheno]]] = pheno
                Intable$Phenotype[Intable$Phenotype %in% PhenoOrder[[pheno]]] = pheno
            }
            pheno_vector = unique(csd$Phenotype)
        } else{
            
            
            PhenoOrder = unique(csd$Phenotype)
            pheno_vector = PhenoOrder
            names(PhenoOrder) = PhenoOrder
            ColsOrder = ColsOrder[names(ColsOrder) %in% PhenoOrder]
            
            
            check_elsestate = TRUE
        }
        colors_phenotype = ColsOrder
        
    }

    if (TRUE %in% plotter){
        output_dir <- file.path(fig.prefix, sample_name)
        if (!dir.exists(output_dir)){
            dir.create(output_dir, recursive = T)
            print(paste("Directory created for", samplename, "with output directory", output_dir))
        } else {
            print(paste("Directory for ", samplename, "already exists with output directory", output_dir, "! Figures were overwritten."))
        }
    }
    
    if (is.null(r_vec)){
        dim_scale = min(max(csd[[XposCol]], max(csd[[YposCol]])))
        r_vec = dim_scale*c(0.1,0.9)
        # r_vec currently does not work well when given NULL due to several automated settings for the grid of r in the settings of each option
        # ASSUMPTION
    }
    
    
    options_all = list("G","F", "J","Gdot", "Jdot", "K", "L", "pcf", "Kdot", "Ldot")
    
    
    if (is.null(options)){
        options = options_all
    } else if (all(options %in% options_all)){
        options = options
    } else{
        stop("one or more spatial statistics in parameter 'options' are not correctly defined")
    }
    
    
    
    # generate pairwise distance matrix for csd for use in getMAD
    
    pairwise_distance = distance_matrix(csd)
    
    
    
    # normal statistics: MAD and MED
    
    # Create Intable with nearest distances for each phenotype, here after the substitution of "" to "Other" and Simplyfying the Phenotypes
    Intable_with_distance = Intable %>%
        do(bind_cols(., find_nearest_distance(.)))
    print(paste0("dimensions of the data with distances is ", dim(Intable_with_distance)[1], " times ", dim(Intable_with_distance)[2]))
    
    output = getMAD(Intable_with_distance, pairwise_distance, pheno_vector)
    MED_min = output[[1]]
    MED = output[[2]]
    MAD_min = output[[3]]
    MAD = output[[4]]
    
    
    # Creation of the Poisson Point Process for plotting overall figure and quadratcounts figures
    
    
    csd_ppp = ppp(x=csd[[XposCol]], y=csd[[YposCol]], 
                  window = owin(c(0, max(csd[[XposCol]])), c(0, max(csd[[YposCol]]))),
                  marks = factor(csd[[PhenoCol]], names(PhenoOrder)))
    unitname(csd_ppp) = c("micron", "microns")
    
    if (plotter[1] == TRUE) {
      
      png(filename = paste0(file.path(output_dir, sample_name),".png"))
      par(mar=rep(0.5, 4))
      #par(mar = c(0,2,0,0)+0.1)
      plot(csd_ppp, cols = colors_phenotype[levels(csd_ppp$marks)], xlab = "", ylab = "", main = "", pch = 20)
      title(paste("Location of cells and their phenotype\n in sample", sample_name), line = -3)
      dev.off()
    }
    
    
    
    
    # normal statistics: Counts and Density
    Area_sample = summary(csd_ppp)$window$area
    output = getDensity(csd, pheno_vector, Area_sample)
    
    counts_sample = output[[1]]
    density_sample = output[[2]]
    
    
    
    
    # Make plots for the quadrat counts of each phenotype and pairwise plots for evaluating the pairwise stats
    
    quadratcount_X2statistic = list()
    quadratcount_X2statistic_normed = list()
    amount_pheno = length(pheno_vector)
    
    for (counter1 in seq_along(pheno_vector)){
      phenotype1 = pheno_vector[counter1]
      sym_matrix_sequence = seq(counter1,seq_along(pheno_vector)[amount_pheno])
      
      for (counter2 in sym_matrix_sequence){
        phenotype2 = pheno_vector[counter2]
        
        splitted = csd_ppp[(marks(csd_ppp) == phenotype1) | (marks(csd_ppp) == phenotype2)] 
        
        if (phenotype1 == phenotype2){
          if (plotter[2] == TRUE){
            # plot quadratcounts and save
            png(filename = paste0(file.path(output_dir, sample_name),"_quadratcounts_", phenotype1, ".png"))
            par(mar=rep(0.5, 4))
            plot(splitted, cols = colors_phenotype[levels(csd_ppp$marks)], xlab = "", ylab = "", main = "",  pch = 20)
            plot(quadratcount(splitted), add = TRUE)
            title(paste("Quadratcounts of", phenotype1, "\n in sample", sample_name), line = -3)
            dev.off()
            # plot singlephenotype and save
            png(filename = paste0(file.path(output_dir, sample_name), '_', phenotype1, ".png"))
            par(mar=rep(0.5, 4))
            plot(splitted, cols = colors_phenotype[levels(csd_ppp$marks)], xlab = "", ylab = "", main = "",  pch = 20)
            title(paste("Location of", phenotype1, "\n in sample", sample_name), line = -3)
            dev.off()
          }
          # normal quadratcounts statistic and normalizing by the total amount of corresponding colour
          
          quadrattest = quadrat.test(splitted) 
          quadratcount_X2statistic[[phenotype1]] = quadrattest$statistic[['X2']]
          quadratcount_X2statistic_normed[[phenotype1]] = quadrattest$statistic[['X2']]/counts_sample[counter1] 
        }
      }
    }
    
    
    
    
    # Replace "+" with "T" and "-" with "F" in all variables for correct functioning of extracting inbuild statistics.
    
    PhenoOrder = lapply(PhenoOrder, function(x) {gsub("+", "T", x, fixed = TRUE)})
    PhenoOrder = lapply(PhenoOrder, function(x) {gsub("-", "F", x, fixed = TRUE)})
    if (isTRUE(check_elsestate)){
      names(PhenoOrder) = PhenoOrder
    }
    
    
    pheno_vector = gsub("+", "T", pheno_vector, fixed = TRUE)
    pheno_vector = gsub("-", "F", pheno_vector, fixed = TRUE)
    
    csd$Phenotype = sapply(csd$Phenotype, function(x) {gsub("+", "T", x, fixed = TRUE)})
    csd$Phenotype = sapply(csd$Phenotype, function(x) {gsub("-", "F", x, fixed = TRUE)})
    
    
    
    csd_ppp = ppp(x=csd[[XposCol]], y=csd[[YposCol]], 
                  window = owin(c(0, max(csd[[XposCol]])), c(0, max(csd[[YposCol]]))),
                  marks = factor(csd[[PhenoCol]], names(PhenoOrder)))
    unitname(csd_ppp) = c("micron", "microns")
    
    values_options = list()
    
    all_types_options_sample_name = list()
    
    r_close_list = list()
    statistic_close_list = list()
    normalized_list = list()
    
    
    for (option in options){
      
      if (option %in% list("K","L","Kdot","Ldot","pcf")){  # (option == "K")|| (option == "L") || (option == "Kdot")|| (option == "Ldot") || (option == "pcf")){
        all_types = alltypes(csd_ppp,fun = paste(option), dataname = sample_name, envelope = envelope_bool, correction = "iso")
      } else {
        all_types = alltypes(csd_ppp,fun = paste(option), dataname = sample_name, envelope = envelope_bool, correction = "km")
      }
      
      all_types_options_sample_name[[option]] = all_types
      
      if (plotter[3] == TRUE){
          
        # if (option %in% list("G","J","K", "L","pcf")){ # (option == "K")|| (option == "L") || (option == "pcf")){
        #     width = 800
        #     height = 700
        #     res = 80
        # } else {
        #     width = 800
        #     height = 900
        #     res = 80
        # }
        # 
        # png(filename = paste0(file.path(output_dir, sample_name),"_statistic_",option,".png"),
        #     width = width, height = height, res = res)
        # plot(all_types, samex = TRUE, samey = TRUE)
        # dev.off()
        
        png(filename = paste0(file.path(output_dir, sample_name),"_statistic_",option,".png"))
        par(mar=rep(0.5, 4))
        plot(all_types, samex = TRUE, samey = TRUE)
        dev.off()
      }
        
      output = interpolate_r(all_types, r_vec, option, envelope_bool)
      
      
      statistic_close_list[[option]] = output[[1]]
      normalized_list[[option]] = output[[2]]
    }
    
    output_all = list()
    output_all[["pairwise_distance"]] = pairwise_distance
    output_all[["counts_sample"]] = counts_sample
    output_all[["Area_sample"]] = Area_sample
    output_all[["density_sample"]] = density_sample
    output_all[["quadratcount_X2statistic"]] = quadratcount_X2statistic
    output_all[["quadratcount_X2statistic_normed"]] = quadratcount_X2statistic_normed
    output_all[["MED_min"]] = MED_min
    output_all[["MED"]] = MED
    output_all[["MAD_min"]] = MAD_min
    output_all[["MAD"]] = MAD
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

# function for receiving the counts of the phenotypes in the sample and their density
getDensity <- function(data, pheno_vector, Area_sample){
    counts_sample = rep(0,length(pheno_vector))
    names(counts_sample) = pheno_vector
    
    density_sample = rep(0,length(pheno_vector))
    names(density_sample) = pheno_vector
    
    for (i in pheno_vector){
        n = dim(data %>% filter(`Phenotype` == i))[1]
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
            } else{
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
                    
                    a = (high2-high1)/(r2-r1) # y=ax+b
                    b = high2 - a*r2
                    high = a*(r_i-r2)+high2
                    
                    lower_bound = all_types[["fns"]][[range]][["lo"]]
                    low1 = lower_bound[left]
                    low2 = lower_bound[right]
                    
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
                    
                    statistic_close_list[[paste("radius", r_i)]][[paste(option, "fns which",range)]] = stat - stat_theo
                    normalized_list[[paste("radius", r_i)]][[paste(option, "fns which",range)]] = normalized
                    
                }
            }
            if ((option == "G")|| (option == "F") || (option == "J")|| (option == "Gdot") || (option == "Jdot")){
                
                
                
                if (max(r_emperic)< r_i){
                    warning("default r interval is to small")
                    next(paste("skip", range, option))
                }
                
                if (envelope_bool == TRUE){
                    
                    higher_bound = all_types[["fns"]][[range]][["hi"]]
                    high1 = higher_bound[left]
                    high2 = higher_bound[right]
                    
                    a = (high2-high1)/(r2-r1) # y=ax+b
                    b = high2 - a*r2
                    high = a*(r_i-r2)+high2
                    
                    lower_bound = all_types[["fns"]][[range]][["lo"]]
                    low1 = lower_bound[left]
                    low2 = lower_bound[right]
                    
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
                    
                    statistic_close_list[[paste("radius", r_i)]][[paste(option, "fns which",range)]] = stat - stat_theo
                    normalized_list[[paste("radius", r_i)]][[paste(option, "fns which",range)]] = normalized
                }
            }
        }
    }
    return(list(statistic_close_list, normalized_list))
}

######### RUN feature extract function #########
feature_extract <- function(outputs, mat_output_dir){
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
    
    
    mat_ripleys = matrix(NA, nrow = length(allfeat_flat), ncol = length(outputs),
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
                    mat_ripleys[df$Ffeatname, name] = df$measure
                    mat_ripleys[paste0('Normalized_', df$Ffeatname), name] = df$Nmeasure
                }
            }
        }
    }
    
    
    counts = c()
    dens = c()
    X2stat = c()
    X2stat_normed = c()
    
    # get phenotypes for counts and densities
    for (out in outputs) {
        counts = union(counts, names(out$counts_sample))
        dens = union(dens, names(out$density_sample))
        X2stat = union(X2stat, names(out$quadratcount_X2statistic))
        X2stat_normed = union(X2stat_normed, names(out$quadratcount_X2statistic_normed))
    }
    
    # create a matrix for counts
    mat_counts = matrix(NA, nrow = length(counts), ncol = length(outputs),
                 dimnames = list(paste0('counts_sample_', sort(counts)), names(outputs)))
    
    # create a matrix for density
    mat_density = matrix(NA, nrow = length(dens), ncol = length(outputs),
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
        data_density = out$density_sample
        data_X2stat = out$quadratcount_X2statistic
        data_X2stat_normed = out$quadratcount_X2statistic_normed
        for (featname in names(data_counts)){
          mat_counts[paste0('counts_sample_',featname),name] = data_counts[[featname]]
          mat_density[paste0('density_sample_',featname),name] = data_density[[featname]]
          mat_X2stat[paste0('X2stat_sample_',featname),name] = data_X2stat[[featname]]
          mat_X2stat_normed[paste0('X2stat_normed_sample_',featname),name] = data_X2stat_normed[[featname]]
        }
    }
    
    
    
    MED_min_pheno = c()
    MED_pheno = c()
    MAD_min_pheno = c()
    MAD_pheno = c()
    
    # get phenotypes for median minimal, median, MAD minimal, MAD.
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
      
      for (row in 1:length(rownames(data_MED))){
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
      
      for (row in 1:length(rownames(data_MAD))){
        featname_from = rownames(data_MAD)[row]
        for (col in row:length(colnames(data_MAD))){
          featname_to = colnames(data_MAD)[col]
          mat_mad[paste0('MAD_', featname_from, '_', featname_to),name] = data_MAD[featname_from, featname_to]
        }
      }
    }
    
    mat = t(rbind(mat_ripleys, mat_counts, mat_density, mat_X2stat, mat_X2stat_normed,
                mat_med_min, mat_med, mat_mad_min, mat_mad))
    
    
    # save matrix to csv for later inspection
    
    write.csv(mat,file = file.path(mat_output_dir, 'feature_matrix.csv'))
    
    
    return(mat)
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
                hist(csd$`Distance to PAX5+PDL1+`[csd$Phenotype == i],main = paste(i,"to PAX5+PDL1+"), ylab = "distance to PAX5+PDL1+")
                boxplot(csd$`Distance to PAX5+PDL1+`[csd$Phenotype == i],main = paste(i,"to PAX5+PDL1+"), ylab = "distance to PAX5+PDL1+")
            }
        }
        
        if ("PAX5+PDL1-" %in% pheno_vector){
            dummy = csd %>%
            filter((`Phenotype` == i))
            statistic_mean_sample_name[i,"Mean shortest distance to PAX5+PDL1-"] = mean(dummy$`Distance to PAX5+PDL1-`)

            if (plotter[1] == TRUE){
                hist(csd$`Distance to PAX5+PDL1-`[csd$Phenotype == i],main = paste(i,"to PAX5+PDL1-"), ylab = "distance to PAX5+PDL1-")
                boxplot(csd$`Distance to PAX5+PDL1-`[csd$Phenotype == i],main = paste(i,"to PAX5+PDL1-"), ylab = "distance to PAX5+PDL1-")
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

# function for calculating area and Maxnorm
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






