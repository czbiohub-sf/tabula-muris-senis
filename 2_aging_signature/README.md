# Aging signiture analysis (Under construction)

We used the **Tabula Muris Senis** dataset to perform a systematic gene-wise study of gene expression changes during aging across different cell types. Stay tuned for the preprint. 

## Data access

The data can be downloaded at [tms_gene_data](https://figshare.com/articles/tms_gene_data/11413869). To run each notebook, please specify the `data_path` variable to be the path of the `tms_gene_data` folder. 

## Data preprocessing 
Code: 
- Preprocessing the FACS data: `filter_FACS.ipynb`
- Preprocessing the droplet data: `filter_droplet.ipynb`

Processed data: 
- The processed FACS data: `tms_gene_data/facs_filtered.h5ad` 
- The processed droplet data: `tms_gene_data/droplet_filtered.h5ad`
- Also bulk data: `tms_gene_data/190304_maca_bulk.h5ad`

## Differential gene expression (DGE) testing

Code: 
- Tissue-cell level DGE testing for the FACS data: `DE_tissue_cell_FACS.ipynb`
- Tissue-cell level DGE testing for the droplet data: `DE_tissue_cell_droplet.ipynb`
- Tissue level DGE testing for the FACS data: `DE_tissue_FACS.ipynb`
- Tissue level DGE testing for the droplet data: `DE_tissue_droplet.ipynb `
- Tissue level DGE testing for the bulk data" `DE_tissue_bulk.ipynb`

DGE testing result: 
- In the folder `tms_gene_data/DE_result`

## Downstream analysis

Code: 
- Tissue-cell level downstream analysis: `downstream_tissue_cell.ipynb`
- Tissue level downstream analysis: `downstream_tissue.ipynb`

Data: 
- The annotation data are in `tms_gene_data/annotation_data`
- The results are in `tms_gene_data/results`
