

library(spatstat)
library(tidyverse)
library(phenoptr)
library(spdep)
library(zoo)

# phenotype can be converted by the "phenotype argument".
# why statistics was used as argument?
do_analyse <- function(Intable, PhenoOrder = NULL, Cols = NULL, phenotype = NULL, plotter = FALSE, 
		       XposCol = 'Cell X Position', YposCol = 'Cell Y Position',
		       PhenoCol = 'Phenotype', ...) {
  
  csd <- Intable[, c(PhenoCol, XposCol, YposCol)]
  colnames(csd) = c('Phenotype', 'Cell X Position',  'Cell Y Position')
  pheno_vector = unique(csd$Phenotype)

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
  

  # print(paste("pheno_vector",pheno_vector))
  # cat("pheno_vector",pheno_vector)
  n=length(pheno_vector)  
  sample_name = str_replace(sample_name,".im3","")
  
  
  # filter first two columns out to convert to ppp as column 1 and 2 are assumed respectively X and Y coordinates and move them in front
  #csd = csd[ ,c(5,6,4,1,2,3,c(7:ncol(csd)))]
  # names(csd[,1:6])
  # [1] "Cell X Position" "Cell Y Position" "Cell ID"         "Path"            "Sample Name"     "Phenotype" 
  
  csd_ppp = ppp(x=csd[['Cell X Position']], y=csd[['Cell Y Position']],
                window = owin(c(0, max(csd[['Cell X Position']])),
                              c(0, max(csd[['Cell Y Position']]))),
                marks = as.factor(csd[['Phenotype']]))
  if (plotter == TRUE){
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
  # colnames(statistics) = c(paste("Area K",pheno_vector_absolut_TF),"Area Kdot", paste("Maxnorm K",pheno_vector_absolut_TF), "Maxnorm Kdot")
  # rownames(statistics) = pheno_vector_absolut_TF  
  values_options = list()
  for (option in options){
    print(option)
    value = alltypes(csd_ppp,fun = paste(option), dataname = sample_name, envelope = FALSE)
    
    values_options[[option]] = value
    if (plotter == TRUE){
      plot(value)
    }
    
    # print(((option == "K") || (option == "Kdot")))
    # if ((option == "K") || (option == "Kdot")){
    #   for (pheno_row in pheno_vector){
    #     # K_self = Kest(csd_ppp[csd_ppp$marks == paste(pheno_row)], correction = "border")
    #     
    #     i = which(pheno_vector %in% pheno_row)
    #     
    #     position = value[[2]][i,i]
    #     
    #     K_self = value[[1]][[position]]
    #     # K_self = Kcs[1+(i-1)*(n+1)][[1]] # [["border"]]
    #     temp = calculate_area_norm(K_self, option)
    #     statistics[paste(pheno_row), paste("Area K",pheno_row)] = temp[1]
    #     statistics[paste(pheno_row), paste("Maxnorm K",pheno_row)] = temp[2]
    #     
    #     if (option == "Kdot"){
    #       # # K_dot = Kdot(csd_ppp,pheno_row, correction = "border")
    #       position = value[[2]][i]
    #       
    #       K_dot = value[[1]][position]
    #       temp = calculate_area_norm(K_dot, option)
    #       statistics[paste(pheno_row), "Area Kdot"] = temp[1]
    #       statistics[paste(pheno_row), "Maxnorm Kdot"] = temp[2]
    #     }
    #     # if (option == "Ldot"){
    #     #   # # L_dot = Ldot(csd_ppp,pheno_row, correction = "border")
    #     #   L_dottemp1 = value[[2]]
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
    #       print(i)
    #       print(j)
    #       print(pheno_col)
    #       position = value[[2]][i,j]  # out of bounds, less phenotypes plotted, maybe because of + en - signs.
    #       print(position)             # gsub("+",".", pheno_vector, fixed = TRUE)
    #       K_cross = value[[1]][[position]] #      Kcs[paste(pheno_row),paste(pheno_col)][["border"]]
    #       
    #       temp = calculate_area_norm(K_cross, option)
    #       statistics[paste(pheno_row), paste("Area K",pheno_col)] = temp[1]
    #       statistics[paste(pheno_row), paste("Maxnorm K",pheno_col)] = temp[2]
    #     }
    #   }
    }
  
  return(c(statistics, values_options))
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
