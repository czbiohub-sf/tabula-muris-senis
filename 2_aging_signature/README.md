# Aging signiture analysis

We used the **Tabula Muris Senis** dataset to perform a systematic gene-wise study of gene expression changes during aging across different cell types. Stay tuned for the preprint. 

- The current document is for the version `revision1`
- The code for the initial submission is in `./archive_initial_submission`

## Data access

The data can be downloaded at [tms_gene_data](https://figshare.com/account/projects/64982/articles/12827615) (the name in the link is `tms_gene_data_rv1`). To run each notebook, please specify the `data_path` variable to be the path of the `tms_gene_data` folder. 

## Datasets
Processed data: 
- The TMS FACS data: `tms_gene_data/tabula-muris-senis-facs-official-raw-obj.h5ad` 
- The TMS mutation data wit ERCC counts: `tms_gene_data/adata_with_ercc_gatk_all_data_with_metadata.h5ad` 
- The TMS droplet data: `tms_gene_data/tabula-muris-senis-droplet-official-raw-obj.h5ad`
- The bulk data [Schaum et al. Nature 2020]: `tms_gene_data/190304_maca_bulk.h5ad`
- The Kimmel et al. data [Kimmel et al. Genome Research 2019]: `tms_gene_data/Kimmel_GR_2019_data.zip` 
- The Kowalczyk et al. data [Kowalczyk et al. Genome Research 2015]: `tms_gene_data/Kowalczyk_GR_2015.zip` 

## Differential gene expression (DGE) testing

Code: 
- Main R script: `tms_gene_data/job.DGE_analysis/DGE_analysis.R`
- DEG analysis with various data and various configurations: `tms_gene_data/job.DGE_analysis/DGE_analysis*.sh`
- Processing script to get the `.rds` data from the `.h5ad` data: `tms_gene_data/job.DGE_analysis/get_data_info.ipynb`

Data:
- Per-tissue `.rds` format TMS data: `tms_gene_data/rds_by_tissue.1e4.zip`
- All DGE results: `tms_gene_data/DGE_result.zip`

## Downstream analysis

Code: 
- Summarizing the DGE results: `tms_gene_data/job.downstream/downstream_tissue_cell.dge_summary.ipynb`
- Partitioning the genes into global aging genes and category-specific aging genes: `tms_gene_data/job.downstream/downstream_tissue_cell.gene_partition.ipynb`
- Gene sets and pathway analysis: `tms_gene_data/job.downstream/downstream_tissue_cell.geneset.ipynb`
- GAG score analysis: `tms_gene_data/job.downstream/downstream_tissue_cell.aging_score.ipynb`

Data: 
- The annotation data: `tms_gene_data/annotation_data.zip`
- All results: `tms_gene_data/result_v1.zip`. 
- Specifically, the gene partition statistics may be of interests to researchers `tms_gene_data/result_v1/tms_gene_table/gene_stats_*.gz`
