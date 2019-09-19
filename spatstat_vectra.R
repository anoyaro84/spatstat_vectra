

library(spatstat)
library(tidyverse)
library(phenoptr)
library(spdep)
library(zoo)

# phenotype can be converted by the "phenotype argument".
# why statistics was used as argument?
do_analyse <- function(Intable, PhenoOrder = NULL, Cols = NULL, phenotype = NULL, plotter = FALSE, 
                       XposCol = 'Cell X Position', YposCol = 'Cell Y Position',
                       PhenoCol = 'Phenotype', sample_name = NULL ...) {
  csd <- Intable[, c(PhenoCol, XposCol, YposCol)]
  colnames(csd) = c('Phenotype', 'Cell X Position',  'Cell Y Position')
  pheno_vector = unique(csd$Phenotype)

  if (is.null(sample_name)) {
	  sample_name = 'Input sample'
  }
  
  if (is.null(PhenoOrder)) {
    PhenoOrder = pheno_vector # if no order is set, just take the order from data
  }
  
  pheno_vector = pheno_vector[order(match(pheno_vector,PhenoOrder))]
  
  # Replace phenotype label by definition
  if (!is.null(phenotype)) {
    for(pheno in names(phenotype)) {
      csd$Phenotype[[pheno]] = phenotype[[pheno]]
    }
  }
  # Replace "+" with "T" and "-" with "F" in the vector and csd
  pheno_vector = str_replace_all(pheno_vector,"[+]","T")
  pheno_vector = str_replace_all(pheno_vector,"[-]","F")
  
  csd$Phenotype = str_replace_all(csd$Phenotype,"[+]","T")
  csd$Phenotype = str_replace_all(csd$Phenotype,"[-]","F")
  
  n=length(pheno_vector)  
  #sample_name = str_replace(sample_name,".im3","")
  
  
  # filter first two columns out to convert to ppp as column 1 and 2 are assumed respectively X and Y coordinates and move them in front
  #csd = csd[ ,c(5,6,4,1,2,3,c(7:ncol(csd)))]
  # names(csd[,1:6])
  # [1] "Cell X Position" "Cell Y Position" "Cell ID"         "Path"            "Sample Name"     "Phenotype" 
  
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
  # options = c("pcf", "G", "Jcross", "Kcross", "Lcross", "Gdot", "Jdot", "Kdot", "Ldot")
  
  
  dim_n = length(PhenoOrder)
  dim_m = 2*dim_n + 2
  print(dim_n)
  
  
  statistics = as.data.frame(matrix(rep(0,dim_n*dim_m),dim_n, dim_m))
  colnames(statistics) = c(paste("Area K", PhenoOrder),"Area Kdot", paste("Maxnorm K", PhenoOrder), "Maxnorm Kdot")
  rownames(statistics) = PhenoOrder
  
  values_options = list()
  
  all_types_options_sample_name = list()
  # print("loop 1: option")
  
  for (option in options){
    # print(option)
    all_types = alltypes(csd_ppp,fun = paste(option), dataname = sample_name, envelope = FALSE)
    # print(length(all_types$fns))
    
    all_types_options_sample_name[[option]] = all_types
    if (plotter[2] == TRUE){
      plot(all_types)
    }
    
    output = interpolate_r(all_types, r_vec, option)
    
    r_close_list = output[1]
    statistic_close_list = output[2]
    
    # print("einde loop 3")
    # print("r_close_list")
    # print(r_close_list)
    # print("statistic_close_list")
    # print(statistic_close_list)
    
    
    
    # print(((option == "K") || (option == "Kdot")))
    # if ((option == "K") || (option == "Kdot")){
    #   for (pheno_row in pheno_vector){
    #     # K_self = Kest(csd_ppp[csd_ppp$marks == paste(pheno_row)], correction = "border")
    #     
    #     i = which(pheno_vector %in% pheno_row)
    #     
    #     position = all_types[[2]][i,i]
    #     
    #     K_self = all_types[[1]][[position]]
    #     # K_self = Kcs[1+(i-1)*(n+1)][[1]] # [["border"]]
    #     temp = calculate_area_norm(K_self, option)
    #     statistics[paste(pheno_row), paste("Area K",pheno_row)] = temp[1]
    #     statistics[paste(pheno_row), paste("Maxnorm K",pheno_row)] = temp[2]
    #     
    #     if (option == "Kdot"){
    #       # # K_dot = Kdot(csd_ppp,pheno_row, correction = "border")
    #       position = all_types[[2]][i]
    #       
    #       K_dot = all_types[[1]][position]
    #       temp = calculate_area_norm(K_dot, option)
    #       statistics[paste(pheno_row), "Area Kdot"] = temp[1]
    #       statistics[paste(pheno_row), "Maxnorm Kdot"] = temp[2]
    #     }
    #     # if (option == "Ldot"){
    #     #   # # L_dot = Ldot(csd_ppp,pheno_row, correction = "border")
    #     #   L_dottemp1 = all_types[[2]]
    #     #   L_dot = L_dottemp1[[paste(pheno_row)]]# [["border"]]
    #     #   temp = calculate_area_norm(L_dot)
    #     #   statistics[paste(pheno_row), "Area Ldot"] = temp[1] # no Ldot in statistics yet
    #     #   statistics[paste(pheno_row), "Maxnorm Ldot"] = temp[2]
    #     # }
    #     pheno_dummy = pheno_vector[pheno_row != pheno_vector]
    #     for (pheno_col in pheno_vector[pheno_row != pheno_vector]){
    #       # K_cross = Kcross(csd_ppp,pheno_row,pheno_col, correction = "border")
    #       
    #       j = which(pheno_vector %in% pheno_col)
    #       position = all_types[[2]][i,j]  # out of bounds, less phenotypes plotted, maybe because of + en - signs.
    #       K_cross = all_types[[1]][[position]] #      Kcs[paste(pheno_row),paste(pheno_col)][["border"]]
    #       
    #       temp = calculate_area_norm(K_cross, option)
    #       statistics[paste(pheno_row), paste("Area K",pheno_col)] = temp[1]
    #       statistics[paste(pheno_row), paste("Maxnorm K",pheno_col)] = temp[2]
    #     }
    #   }
  }
  
  return(output)
}
# function for receiving the counts of the phenotypes in the sample and their density
getDensity <- function(data, pheno_vector,phenoOrder, Area_sample){
  counts_sample = rep(0,length(phenoOrder))
  names(counts_sample) = phenoOrder
  
  density_sample = rep(0,length(phenoOrder))
  names(density_sample) = phenoOrder
  
  for (i in pheno_vector){
    n = dim(data %>% filter(`Phenotype` == i))[1]
    counts_sample[[i]] = n
    density_sample[[i]] = n / Area_sample
  }
  
  return(c(counts_sample, density_sample))
}

# function interpolate statistic for given r
interpolate_r <- function(all_types, r_vec, option){
  for (r_i in r_vec){
    r_close_list[[option]][[paste("radius", r_i)]] = list()
    statistic_close_list[[option]][[paste("radius", r_i)]] = list()
    
    for (range in (1:length(all_types$fns))){
      # print("inside loops, test r_emperic")
      # print(all_types_options[[option]][["fns"]][[range]][["r"]])
      r_emperic = all_types[["fns"]][[range]][["r"]]
      
      dif = abs(r_emperic - r_i)
      
      index1 = match(min(dif),dif)
      
      index2 = index1 + 1
      r1 = r_emperic[index1]
      r2 = r_emperic[index2]
      # print(c(r1,r2))
      
      
      r_close_list[[option]][[paste("radius", r_i)]][[paste(option, "fns which",range)]] = r_i
      
      
      if ((option == "K")|| (option == "L") || (option == "Kdot")|| (option == "Ldot")){
        stat1 = all_types[["fns"]][[range]][["border"]][index1]
        stat2 = all_types[["fns"]][[range]][["border"]][index2]
        
        a = (stat2-stat1)/(r2-r1)
        b = stat2 - a*r2
        stat = a*(r_i-r2)+stat2
        statistic_close_list[[option]][[paste("radius", r_i)]][[option]][[paste(option, "fns which",range)]] = stat
      } else if (option == "pcf"){
        stat1 = all_types[["fns"]][[range]][["trans"]][index1]
        stat2 = all_types[["fns"]][[range]][["trans"]][index2]
        
        a = (stat2-stat1)/(r2-r1)
        b = stat2 - a*r2
        stat = a*(r_i-r2)+stat2
        
        statistic_close_list[[option]][[paste("radius", r_i)]][[option]][[paste(option, "fns which",range)]] = stat
      } else{
        stat1 = all_types[["fns"]][[range]][["rs"]][index1]
        stat2 = all_types[["fns"]][[range]][["rs"]][index2]
        
        a = (stat2-stat1)/(r2-r1)
        b = stat2 - a*r2
        stat = a*(r_i-r2)+stat2
        statistic_close_list[[option]][[paste("radius", r_i)]][[option]][[paste(option, "fns which",range)]] = stat
      }
    }
  }
  return(c(r_close_list, statistic_close_list))
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
