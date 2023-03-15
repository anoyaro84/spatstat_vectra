

# parse input arguments
if (exists('snakemake')) { # if the script is used by snakemake
    infile = snakemake@input[["infile"]]
    outRDS = snakemake@output[["RDS"]]
    r_vec = as.numeric(strsplit(snakemake@params[["r_vec"]], split=',')[[1]])
    #pheno_vector_absolut = strsplit(snakemake@params[["phenotype"]], split=',')[[1]]
    pheno_vector_absolut = snakemake@params[["phenotype"]]
    colors_absolut = strsplit(snakemake@params[["colors"]], split=',')[[1]]
    #colors_absolut = snakemake@params[["colors"]]
    path_script = 'scripts/spatstat_vectra'
    path_figure = snakemake@params[["fig_prefix"]]
    plotOnly <- snakemake@params[["plotOnly"]]
    ref_ctype = strsplit(snakemake@params[["ref_ctype"]], split=',')[[1]]
}else {
    args = commandArgs(trailingOnly=TRUE)
    infile = args[1]
    outRDS = args[2]
    r_vec = as.numeric(strsplit(args[3], split=',')[[1]])
    pheno_vector_absolut = strsplit(args[4], split=',')[[1]]
    color_absolut = strsplit(args[5], split=',')[[1]]
    ref_ctype = strsplit(args[6], split=',')[[1]]
    path_script = 'scripts/spatstat_vectra'
    path_figure = args[8]


    
    #infile <- 'data/vectra/processed/Lympho/TVU07-18235 I-C_cell_seg_data.txt'
    #r_vec <- as.numeric(strsplit('5,10,25,50,75,100,150,200,250,500', split=',')[[1]])
    #pheno_vector_absolut <- c('T cell Other','T cell Other (memory)','Cytotoxic T cell','Cytotoxic T cell (memory)','T helper cell','T helper cell (memory)','Regulatory T cell','Regulatory T cell (memory)','B cell','Tumor cell','DAPI','Unknown cell')
    #colors_absolut <- c('#800000','#e6194B','#9A6324','#f58231','#ffe119','#fffac8','#3cb44b','#aaffc3','#000075','#000000','#911eb4','#a9a9a9')
    #ref_ctype <- 'Tumor cell'
    
}

if (is.character(pheno_vector_absolut)) {
    pheno_vector_absolut = strsplit(pheno_vector_absolut, split=',')[[1]]
}

print(r_vec)
print(pheno_vector_absolut)
print(colors_absolut)
print(file.path(path_script, 'spatstat_vectra.R'))

names(pheno_vector_absolut) <- pheno_vector_absolut
# load packages
source(file.path(path_script, 'spatstat_vectra.R'))


# processing segmentation file
print(paste('processing', infile))

# extract sample name
samplename = gsub('_cell_seg_data.txt', '', tail(strsplit(infile, '/')[[1]],1))

# set seed for reproducing
set.seed(0)

# perform analysis on segmentation file
output = do_analyse(seg_path = infile, PhenoOrder=pheno_vector_absolut, ColsOrder=colors_absolut, 
			XposCol = 'Cell X Position', YposCol = 'Cell Y Position', PhenoCol = 'Phenotype',
			sample_name = samplename, plotter = c(TRUE, TRUE,TRUE), fig.prefix = path_figure,
			r_vec = r_vec, spatstat_statistics = 'ALL', reference=ref_ctype,plotOnly=plotOnly)

print(output)
print(length(output))
# save outcome
saveRDS(output, outRDS)

