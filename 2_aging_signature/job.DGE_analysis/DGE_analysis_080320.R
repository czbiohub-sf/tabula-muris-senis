#!/usr/bin/env Rscript

start_time =proc.time()

library(MAST)

# Load parameters 
args = commandArgs(trailingOnly=TRUE)
output_folder = args[1]
data_name = args[2]
str_n_gene = args[3]
analyte = args[4]
tissue = unlist(strsplit(analyte, "[.]"))[1]
celltype = unlist(strsplit(analyte, "[.]"))[2]
celltype = gsub('_', ' ', celltype)

print(paste0('output_folder: ', output_folder))
print(paste0('str_n_gene: ', str_n_gene))
print(paste0('tissue: ', tissue))
print(paste0('celltype: ', celltype))

# Load data 
adata_temp = readRDS(paste0('/n/groups/price/martin/tms_gene_data/rds_by_tissue.1e4/',
                            data_name, '.normalized.', tissue, '.rds'))

print('Before filtering')
print(dim(adata_temp))

if (is.na(celltype)==FALSE){
    print('tc analysis')
    ind_select = (adata_temp$cell_ontology_class==celltype)
    adata_temp = adata_temp[,ind_select]
} else {
    print('tissue analysis')
}

print('After filtering')
print(dim(adata_temp))

# Prepare sca object
sca <- SceToSingleCellAssay(adata_temp, class = "SingleCellAssay")
colData(sca)$n_genes = scale(colData(sca)$n_genes) # n_gene (CDR)
sca_filt = sca[rowSums(assay(sca)) != 0, ]

# Set flags 
if (length(unique(sca_filt$sex))>1){
    flag_sex=TRUE
} else {
    flag_sex=FALSE
}
print(paste0('flag_sex: ', flag_sex))
if (str_n_gene=='T'){
    flag_n_gene=TRUE
} else {
    flag_n_gene=FALSE
}
print(paste0('flag_n_gene: ', flag_n_gene))

# DGE testing 
covariate = ''
if (flag_sex==TRUE){
    covariate = paste0(covariate, " + sex")
}
if (flag_n_gene==TRUE){
    covariate = paste0(covariate, " + n_genes")
}

print(paste0('covariate: ', covariate))

zlmCond <- zlm(formula = as.formula(paste0("~age_num", covariate)), sca=sca_filt)

summaryCond <- summary(zlmCond, doLRT="age_num")
if (flag_sex==TRUE){
    summaryCond_sex <- summary(zlmCond, doLRT='sexmale')
}

# Summarize results 
summaryDt <- summaryCond$datatable
dt1 = summaryDt[contrast=="age_num" & component=="H", .(primerid, `Pr(>Chisq)`)]
dt2 = summaryDt[contrast=="age_num" & component=="logFC", .(primerid, coef, z)]
de_res = merge(dt1, dt2, by="primerid")
colnames(de_res) <- c("gene", "age.H_p", "age.logFC", 'age.logFC_z')
de_res$age.H_fdr <- p.adjust(de_res$age.H_p, "fdr")

if (flag_sex==TRUE){
    summaryDt <- summaryCond_sex$datatable
    dt_sex1 = summaryDt[contrast=="sexmale" & component=="H", .(primerid, `Pr(>Chisq)`)]
    dt_sex2 = summaryDt[contrast=="sexmale" & component=="logFC", .(primerid, coef, z)]
    de_res_sex = merge(dt_sex1, dt_sex2, by="primerid")
    colnames(de_res_sex) <- c("gene", "sexmale.H_p", "sexmale.logFC", 'sexmale.logFC_z')
    de_res_sex$sexmale.H_fdr <- p.adjust(de_res_sex$sexmale.H_p, "fdr")
    de_res = merge(de_res, de_res_sex, by="gene")
}

# # Summarize results 
# summaryDt <- summaryCond$datatable
# dt1 = summaryDt[contrast=="age_num" & component=="H", .(primerid, `Pr(>Chisq)`)]
# dt2 = summaryDt[contrast=="age_num" & component=="logFC", .(primerid, coef)]
# de_res = merge(dt1, dt2, by="primerid")
# colnames(de_res) <- c("gene", "age.raw_p", "age.coef")
# de_res$age.bh_p <- p.adjust(de_res$age.raw_p, "fdr")
# 
# if (flag_sex==TRUE){
#     summaryDt <- summaryCond_sex$datatable
#     dt_sex1 = summaryDt[contrast=="sexmale" & component=="H", .(primerid, `Pr(>Chisq)`)]
#     dt_sex2 = summaryDt[contrast=="sexmale" & component=="logFC", .(primerid, coef)]
#     de_res_sex = merge(dt_sex1, dt_sex2, by="primerid")
#     colnames(de_res_sex) <- c("gene", "sexmale.raw_p", "sexmale.coef")
#     de_res_sex$sexmale.bh_p <- p.adjust(de_res_sex$sexmale.raw_p, "fdr")
#     de_res = merge(de_res, de_res_sex, by="gene")
# }

# Write results
write.csv(de_res,paste0(output_folder, '/',analyte,'.csv'), row.names=FALSE)

print('Finished')
print(proc.time() - start_time)

