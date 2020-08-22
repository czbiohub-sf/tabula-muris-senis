# Aging signiture analysis

We used the **Tabula Muris Senis** dataset to perform a systematic gene-wise study of gene expression changes during aging across different cell types. 

Please see more information in the biorxiv preprint Zhang et al. "Mouse Aging Cell Atlas Analysis Reveals Global and Cell Type Specific Aging Signatures" [link](https://www.biorxiv.org/content/10.1101/2019.12.23.887604v2)

- The current document is for the version `revision1`
- The code for the initial submission is in `./archive_initial_submission`

## Data access

The data can be downloaded at [tms_gene_data](https://figshare.com/account/projects/64982/articles/12827615) (the name in the link is `tms_gene_data_rv1`). To run each notebook, please specify the `DATA_PATH` variable to be the path of the `tms_gene_data` folder. 

## Datasets
Processed data: 
- The TMS FACS data: `tms_gene_data/tabula-muris-senis-facs-official-raw-obj.h5ad` 
- The TMS FACS mutation data with ERCC counts: `tms_gene_data/adata_with_ercc_gatk_all_data_with_metadata.h5ad` 
- The TMS droplet data: `tms_gene_data/tabula-muris-senis-droplet-official-raw-obj.h5ad`
- The bulk data [Schaum et al. Nature 2020]: `tms_gene_data/190304_maca_bulk.h5ad`
- The Kimmel et al. data [Kimmel et al. Genome Research 2019]: `tms_gene_data/Kimmel_GR_2019_data.zip` 
- The Kowalczyk et al. data [Kowalczyk et al. Genome Research 2015]: `tms_gene_data/Kowalczyk_GR_2015.zip` 

## Differential gene expression (DGE) testing

Code: 
- Main R script: `tms_gene_data/job.DGE_analysis/DGE_analysis.R`
- DEG analysis with various data and various configurations: `tms_gene_data/job.DGE_analysis/DGE_analysis*.sh`
- Processing script to get the per-tissue `.rds` data from the `.h5ad` data: `tms_gene_data/job.DGE_analysis/get_data_info.ipynb`

Data:
- Per-tissue `.rds` format TMS data: `tms_gene_data/rds_by_tissue.1e4.zip`
- All DGE results: `tms_gene_data/DGE_result.zip` The columns are 
  1) gene
  2) age.H_p: the p-value from the MAST hurdle model. This was used to compute the FDR in the paper.
  3) age.logFC: log-fold changes (in the unit of logFC per month) from hurdle model components. This was as the age coefficient in the paper.
  4) age.logFC_z: the z-score for age.logFC. This was used to compute the standard error of the age coefficients in the paper (namely, logFC_se=logFC/logFC_z)
  5) age.H_fdr: the FDR corresponding to age.H_p, 
  6) sexmale.H_p, 
  7) sexmale.logFC, 
  8) sexmale.logFC_z, 
  9) sexmale.H_fdr:  

## Downstream analysis

Code: 
- Summarizing the DGE results: `tms_gene_data/job.downstream/downstream_tissue_cell.dge_summary.ipynb`
- Partitioning the genes into global aging genes and category-specific aging genes: `tms_gene_data/job.downstream/downstream_tissue_cell.gene_partition.ipynb`
- Gene sets and pathway analysis: `tms_gene_data/job.downstream/downstream_tissue_cell.geneset.ipynb`
- GAG score analysis: `tms_gene_data/job.downstream/downstream_tissue_cell.aging_score.ipynb`

Data: 
- The annotation data: `tms_gene_data/annotation_data.zip`
- All results: `tms_gene_data/result_v1.zip`. 
- Specifically, the gene partition statistics may be of interests to researchers `tms_gene_data/result_v1/tms_gene_table/gene_stats_*.gz`. The columns are:
  1) gene
  2) prop_sig: proportion of tissue-cell types that the gene is significantly related to aging.
  3) prop_upreg: proportion of tissue-cell types that the gene is up-regulated during aging.
  4) median_fc: median of the log fold change over all tissue-cell types.
  5) median_fdr: median of the FDR (FDR adjusted p-value) over all tissue-cell types.
  6) prop_sig_w: weighted proportion of tissue-cell types that the gene is significantly related to aging. This was used in the paper.
  7) prop_upreg_w: weighted proportion of tissue-cell types that the gene is up-regulated during aging. This was used in the paper.
  8) global: boolean flag indicating if the gene is a GAG. This is consistent with Supp. Table. 3
  9) global.dir: direction of the GAG (if the gene is a GAG). This is consistent with Supp. Table. 3
  10) spec_*.mean: mean of the within-set meta age coefficient. This is consistent with Supp. Table. 3
  11) spec_*.se: se of the within-set meta age coefficient. This is consistent with Supp. Table. 3
  12) spec_*.mean_ref: mean of the outside-set meta age coefficient. This is consistent with Supp. Table. 3
  13) spec_*.se_ref: se of the outside-set meta age coefficient. This is consistent with Supp. Table. 3
  14) spec_*.p_dif: p-value for comparing the within-set meta age coefficient and the outside-set meta age coefficient. This is consistent with Supp. Table. 3
  15) spec_*.fdr_dif: FDR (FDR adjusted p-value) for comparing the within-set meta age coefficient and the outside-set meta age coefficient, with respect to all genes. This is consistent with Supp. Table. 3

## Supp. Tables in the paper
- Supp. Table 1: Summary of tissue and cell types.
- Supp. Table 2: Aging-related genes discovered in each tissue-cell type.
- Supp. Table 3: Genesets identified in the paper.
  1) global_aging: global aging genes.
  2) spec_celltype*: genes specific to a cell type.
  3) spec_func* : genes specific to a functional category.
  4) spec_tissue* : genes specific to a tissue.
