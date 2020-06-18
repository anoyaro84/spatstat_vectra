

# parse input arguments
if (exists('snakemake')) { # if the script is used by snakemake
    infile = snakemake@input[["infile"]]
    outRDS = snakemake@output[["RDS"]]
    r_vec = as.numeric(strsplit(snakemake@params[["r_vec"]], split=',')[[1]])
    #pheno_vector_absolut = strsplit(snakemake@params[["phenotype"]], split=',')[[1]]
    pheno_vector_absolut = snakemake@params[["phenotype"]]
    #colors_absolut = strsplit(snakemake@params[["colors"]], split=',')[[1]]
    colors_absolut = snakemake@params[["colors"]]
    path_script = snakemake@params[["pathscript"]]
    path_figure = snakemake@params[["fig_prefix"]]
}else {
    args = commandArgs(trailingOnly=TRUE)
    infile = args[1]
    outRDS = args[2]
    r_vec = as.numeric(strsplit(args[3], split=',')[[1]])
    pheno_vector_absolut = strsplit(args[4], split=',')[[1]]
    color_absolut = strsplit(args[5], split=',')[[1]]
    path_script = args[6]
    path_figure = args[7]
}


print(r_vec)
print(pheno_vector_absolut)
print(colors_absolut)
print(file.path(path_script, 'spatstat_vectra.R'))

# load packages
source(file.path(path_script, 'spatstat_vectra.R'))


# processing segmentation file
print(paste('processing', infile))

# extract sample name
samplename = gsub('_cell_seg_data.txt', '', tail(strsplit(infile, '/')[[1]],1))

# perform analysis on segmentation file
output = do_analyse(segmentation_path = infile, PhenoOrder=pheno_vector_absolut, ColsOrder=colors_absolut, 
			XposCol = 'Cell X Position', YposCol = 'Cell Y Position', PhenoCol = 'Phenotype',
			sample_name = samplename, plotter = c(TRUE, TRUE,TRUE), fig.prefix = path_figure,
			r_vec = r_vec, spatstat_statistics = 'ALL')


# save outcome
saveRDS(output, outRDS)

