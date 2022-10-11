import glob
from os import path

#PATH_DATA = 'HO105_Component data+segmentation map/'
PATH_DATA = 'new_data/combined/'
PATH_OUT = 'spatstat_RDS/'
PATH_FIG = 'Plots/'
PATH_SCRIPT = 'spatstat_vectra'


Files = glob.glob(path.join(PATH_DATA, '*_data.txt'))
print(Files)
IDs = [a.split('/')[-1].split('_cell')[0] for a in Files]

Radius = '5,10,25,50,75,100,150,200,250,500'

print(IDs)
print(len(IDs))

Tumors = {
    'Bsimple' : {
        'Tumors': ['PAX5+PDL1-', 'PAX5+PDL1+']
    },
    'Bdetail' : {
        'PDL1-Tumors': ['PAX5+PDL1-'],
        'PDL1+Tumors': ['PAX5+PDL1+']
    }
}

Tcells = {
    'Tsimple' : {
        'Tcells': ["CD3+CD8-PD1-", "CD3+CD8-PD1+", "CD3+CD8+PD1-", "CD3+CD8+PD1+"]
    },
    'Tdetail': {
        'CD8+Tcells': ["CD3+CD8+PD1-", "CD3+CD8+PD1+"],
        'CD8-Tcells': ["CD3+CD8-PD1-", "CD3+CD8-PD1+"]
    },
    'TdetailPD1': {
        'CD8+PD1+Tcells': ["CD3+CD8+PD1+"],
        'CD8+PD1-Tcells': ["CD3+CD8+PD1-"],
        'CD8-PD1+Tcells': ["CD3+CD8-PD1+"],
        'CD8-PD1-Tcells': ["CD3+CD8-PD1-"]
    }

}

Macrophage = {
    'Msimple' : {
        'Macrophage': ['CD163+PDL1+', 'CD163+PDL1-']
    },
    'Mdetail' : {
        'PDL1+Macrophage' : ['CD163+PDL1+'],
        'PDL1-Macrophage' : ['CD163+PDL1-']
    }
}

def phenotype(wildcard):
    TumorKey = wildcard.split('_')[0]
    TcellKey = wildcard.split('_')[1]
    MacroKey = wildcard.split('_')[2]

    phenotype = {}
    for d in [Tumors[TumorKey], Tcells[TcellKey], Macrophage[MacroKey]]:
        phenotype.update(d)
    phenotype['Others'] = ['Other', 'Other PDL1+']

    return phenotype

def colors(wildcard):
    TumorKey = wildcard.split('_')[0]
    TcellKey = wildcard.split('_')[1]
    MacroKey = wildcard.split('_')[2]

    colors_all = ["magenta", "brown", "red", "blue", "green", "yellow", "gray", "pink", "orange", "cyan"]

    colors = {}
    ind = 0
    for d in [Tumors[TumorKey], Tcells[TcellKey], Macrophage[MacroKey]]:
        for key in d:
            colors[key] = colors_all[ind]
            ind = ind + 1
    colors['Others'] = colors_all[ind]

    return colors

def refcells(wildcards):
    TumorKey = wildcards.split('_')[0]
    return ','.join(list(Tumors[TumorKey].keys()))

print(phenotype('Bsimple_Tsimple_Msimple'))
print(phenotype('Bsimple_Tdetail_Msimple'))
print(colors('Bsimple_Tdetail_Msimple'))
print(refcells('Bdetail_Tdetail_Msimple'))


rule all:
    input:
        expand(path.join(PATH_OUT, '{sample}.RDS'), sample = IDs),
#        expand(path.join(PATH_OUT, 'simplified', '{sample}.RDS'), sample = IDs),
        #expand(path.join(PATH_OUT, '{tumor}_{tcell}_{macro}', '{sample}.RDS'),
        #    tumor=['Bsimple', 'Bdetail'], tcell=['Tsimple','Tdetail'],
        #    macro=['Msimple', 'Mdetail'], sample = IDs),
        expand(path.join(PATH_OUT, '{tumor}_{tcell}_{macro}', '{sample}.RDS'),
            tumor=['Bsimple'], tcell=['TdetailPD1'],
            macro=['Msimple'], sample = IDs),


rule spatstat_general:
    input:
        infile=path.join(PATH_DATA, '{sample}_cell_seg_data.txt')
    output:
        RDS=path.join(PATH_OUT, '{classification}', '{sample}.RDS')
    params:
        #r_vec="2,3,5,10,20,30,40,50",
        r_vec=Radius,
        fig_prefix='{classification}_Figures/',
        phenotype = lambda wildcards: phenotype(wildcards.classification),
        colors= lambda wildcards: colors(wildcards.classification),
    	pathscript = PATH_SCRIPT,
        ref_ctype = lambda wildcards: refcells(wildcards.classification)
    script:
        path.join(PATH_SCRIPT, 'spatstat_commandline.R')



rule spatstat_simple:
    input:
        infile=path.join(PATH_DATA, '{sample}_cell_seg_data.txt')
    output:
        RDS=path.join(PATH_OUT, 'simplified', '{sample}.RDS')
    params:
        #r_vec="2,3,5,10,20,30,40,50",
        r_vec=Radius,
        fig_prefix='SimplfiedFigures/',
        phenotype= {
            "Macrophage": ["CD163+PDL1-","CD163+PDL1+"],
            "Tcells": ["CD3+CD8-PD1-", "CD3+CD8-PD1+", "CD3+CD8+PD1-", "CD3+CD8+PD1+"],
            "Tumors": ["PAX5+PDL1-", "PAX5+PDL1+"],
            "Others": ["Other","Other PDL1+"]
        },
        colors={            
            "Macrophage": "magenta",
            "Tcells": "red",
            "Tumors": "orange",
            "Others": "gray"
        },
    	pathscript = PATH_SCRIPT,
        ref_ctype = 'Tumors'
    script:
        path.join(PATH_SCRIPT, 'spatstat_commandline.R')

rule spatstat:
    input:
        infile=path.join(PATH_DATA, '{sample}_cell_seg_data.txt')
    output:
        RDS=path.join(PATH_OUT, '{sample}.RDS')
    params:
        #r_vec="2,3,5,10,20,30,40,50",
        r_vec=Radius,
        fig_prefix='Figures/',
        phenotype=','.join(["CD163+PDL1-","CD163+PDL1+","CD3+CD8-PD1-", "CD3+CD8-PD1+", "CD3+CD8+PD1-", "CD3+CD8+PD1+","Other","Other PDL1+","PAX5+PDL1-", "PAX5+PDL1+"]),
        colors=','.join(["magenta", "brown", "red", "blue", "green", "yellow", "gray", "pink", "orange", "cyan"]),
        pathscript = PATH_SCRIPT,
        ref_ctype = ','.join(["PAX5+PDL1-", "PAX5+PDL1+"])
    script:
        path.join(PATH_SCRIPT, 'spatstat_commandline.R')
