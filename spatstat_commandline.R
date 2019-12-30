

# parse input arguments
if (exists(snakemake)) { # if the script is used by snakemake
    infile = snakemake@input[["infile"]]
    outRDS = snakemake@output[["RDS"]]
    r_vec = as.numeric(snakemake@params[["r_vec"]])
    pheno_vector_absolut = snakemake@params[["phenotype"]]
    colors_absolut = snakemake@params[["colors"]]
}else {
    args = commandArgs(trailingOnly=TRUE)
    infile = args[1]
    outRDS = args[2]
    r_vec = as.numeric(strsplit(args[3], split=',')[[1]])
    pheno_vector_absolut = strsplit(args[4], split=',')[[1]]
    color_absolut = strsplit(args[5], split=',')[[1]]
}


# load packages
library(reshape2)
source('spatstat_vectra.R')


# load data
print(paste('procesing', infile))

data = purrr::map_df(infile, read_cell_seg_data, pixels_per_micron = getOption("phenoptr.pixels.per.micron"))
data_with_distance = data %>%
	  do(bind_cols(., find_nearest_distance(.)))

# extract sample name
samplename = gsub('_cell_seg_data.txt', '', tail(strsplit(infile, '/')[[1]],1))

# perform analysis
output = do_analyse(data_with_distance, pheno_vector_absolut, colors_absolut, 
			NULL, plotter = c(TRUE, TRUE,TRUE),XposCol = 'Cell X Position', envelope=TRUE, 
			fig.prefix = './', fig.width = 720, fig.height = 720,
			YposCol = 'Cell Y Position',PhenoCol = 'Phenotype', samplename, r_vec = r_vec)


# save outcome
saveRDS(output, outRDS)
