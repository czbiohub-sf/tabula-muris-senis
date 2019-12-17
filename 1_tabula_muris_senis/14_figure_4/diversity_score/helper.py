import scanpy.api as sc
import pandas as pd
import numpy as np
from anndata import read_h5ad
from anndata import AnnData
from itertools import product

def load_normalized_data(file_path, log1p=True):
    """load normalized data
    1. Load filtered data for both FACS and droplet
    2. Size factor normalization to counts per 10 thousand
    3. log(x+1) transform
    4. Combine the data 

    Args:
        file_path (str): file path.

    Returns:
        adata_combine (AnnData): Combined data for FACS and droplet
    """
    # Load filtered data
    # adata_facs = read_h5ad(f'{file_path}/facs_filtered.h5ad')
    adata_facs = read_h5ad(f'{file_path}/facs_filtered_reannotated-except-for-marrow-lung-kidney.h5ad')
    adata_droplet = read_h5ad(f'{file_path}/droplet_filtered.h5ad')
    # Size factor normalization
    sc.pp.normalize_per_cell(adata_facs, counts_per_cell_after=1e4)
    sc.pp.normalize_per_cell(adata_droplet, counts_per_cell_after=1e4)
    # log(x+1) transform
    if log1p:
        sc.pp.log1p(adata_facs)
        sc.pp.log1p(adata_droplet)
    # Combine the data 
    ind_select = adata_facs.obs['age'].isin(['3m', '18m', '24m'])
    adata_facs = adata_facs[ind_select,]
    adata_combine = AnnData.concatenate(adata_facs, adata_droplet, batch_key='b_method',
                                        batch_categories = ['facs','droplet'])
    return adata_combine

def load_normalized_data_full(file_path, log1p=True):
    """load normalized data (full data, including 18m FACS data)
    1. Load filtered data for both FACS and droplet
    2. Size factor normalization to counts per 10 thousand
    3. log(x+1) transform
    4. Combine the data 

    Args:
        file_path (str): file path.
        log1p (bool): If to perform log1o transform

    Returns:
        adata_combine (AnnData): Combined data for FACS and droplet
    """
    # Load filtered data
    adata_facs = read_h5ad(f'{file_path}/facs_filtered_full.h5ad')
    adata_droplet = read_h5ad(f'{file_path}/droplet_filtered.h5ad')
    # Size factor normalization
    sc.pp.normalize_per_cell(adata_facs, counts_per_cell_after=1e4)
    sc.pp.normalize_per_cell(adata_droplet, counts_per_cell_after=1e4)
    # log(x+1) transform
    if log1p:
        sc.pp.log1p(adata_facs)
        sc.pp.log1p(adata_droplet)
    # Combine the data 
    ind_select = adata_facs.obs['age'].isin(['3m', '18m', '24m'])
    adata_facs = adata_facs[ind_select,]
    adata_combine = AnnData.concatenate(adata_facs, adata_droplet, batch_key='b_method',
                                        batch_categories = ['facs','droplet'])
    return adata_combine

def load_normalized_data_bulk(file_path, log1p=True):
    """load normalized bulk data
    1. Load filtered data for bulk RNA-Seq
    2. Size factor normalization to counts per 10 thousand
    3. log(x+1) transform

    Args:
        file_path (str): file path.
        log1p (bool): If to perform log1o transform

    Returns:
        adata_combine (AnnData): Combined data for FACS and droplet
    """
    # Load filtered data
    adata_bulk = read_h5ad(file_path+'/Bulk_Data/190304_maca_bulk.h5ad')
    # QC filtering
    sc.pp.filter_genes(adata_bulk, min_cells=5)
    sc.pp.filter_cells(adata_bulk, min_genes=500)
    adata_bulk = adata_bulk[adata_bulk.obs.Organ!='nan', :]
    # Additional annotations
    adata_bulk.obs['n_counts'] = np.sum(adata_bulk.X, axis=1)
    mito_genes = adata_bulk.var_names.str.startswith('mt-')
    adata_bulk.obs['percent_mito'] = np.sum(adata_bulk[:, mito_genes].X, axis=1) / \
                                        np.sum(adata_bulk.X, axis=1)
    adata_bulk.obs['tissue'] = adata_bulk.obs['Organ']
    adata_bulk.obs['age'] = ['%dm'%x for x in adata_bulk.obs['Age']]
    adata_bulk.obs['age_num'] = [int(x) for x in adata_bulk.obs['Age']]
    adata_bulk.obs['sex'] = ['male' if x=='m' else 'female' for x in adata_bulk.obs['Sex']]
    # Size factor normalization
    sc.pp.normalize_per_cell(adata_bulk, counts_per_cell_after=1e4)
    # log(x+1) transform
    if log1p:
        sc.pp.log1p(adata_bulk)
    # # Combine the data 
    # ind_select = adata_facs.obs['age'].isin(['3m', '18m', '24m'])
    # adata_facs = adata_facs[ind_select,]
    # adata_combine = AnnData.concatenate(adata_facs, adata_droplet, batch_key='b_method',
    #                                     batch_categories = ['facs','droplet'])
    return adata_bulk

# call MAST for DE analysis: numerical
def call_MAST_age():
    """Return the cmd to call MAST to test for age 
    1. Test for age (numerical)
    2. Control for sex

    Args:

    Returns:
        R_cmd (str): R cmd to call MAST
    """
    R_cmd_list = ['sca <- SceToSingleCellAssay(adata_temp, class = "SingleCellAssay")',
                  'sca_filt = sca[rowSums(assay(sca)) != 0, ]',
                  'zlmCond <- zlm(formula = as.formula(paste0("~condition", covariate)),',
                  'sca=sca_filt)',
                  'summaryCond <- summary(zlmCond, doLRT="condition")',
                  'summaryDt <- summaryCond$datatable',
                  'de_res <- merge(summaryDt[contrast=="condition" & component=="H",',
                  '.(primerid, `Pr(>Chisq)`)], #P-vals',
                  'summaryDt[contrast=="condition" & component=="logFC",',
                  '.(primerid, coef)],',
                  'by="primerid") #logFC coefficients',
                  'de_res <- de_res[,FDR:=p.adjust(`Pr(>Chisq)`, "fdr")]']
    R_cmd = '\n'.join(R_cmd_list)
    return R_cmd

# call MAST for DE analysis for sex
def call_MAST_sex():
    """Return the cmd to call MAST to test for sexual dimorphism 
    1. Male (case) v.s. female (control)
    2. age_num (numerical) and n_genes (numerical) as covariates

    Args:

    Returns:
        R_cmd (str): R cmd to call MAST
    """
    R_cmd_list = ['sca <- SceToSingleCellAssay(adata_temp, class = "SingleCellAssay")',
                  'colData(sca)$n_genes = scale(colData(sca)$n_genes)',
                  'sca_filt = sca[rowSums(assay(sca)) != 0, ]',
                  'cond<-factor(colData(sca_filt)$sex)',
                  'cond<-relevel(cond, "female")',
                  'colData(sca_filt)$sex <- cond',
                  'zlmCond <- zlm(formula = ~sex + n_genes + age_num, sca=sca_filt)',
                  'summaryCond <- summary(zlmCond, doLRT="sexmale")',
                  'summaryDt <- summaryCond$datatable',
                  'de_res <- merge(summaryDt[contrast=="sexmale" & component=="H",',
                  '.(primerid, `Pr(>Chisq)`)], #P-vals',
                  'summaryDt[contrast=="sexmale" & component=="logFC",',
                  '.(primerid, coef)],',
                  'by="primerid") #logFC coefficients',
                  'de_res <- de_res[,FDR:=p.adjust(`Pr(>Chisq)`, "fdr")]']
    R_cmd = '\n'.join(R_cmd_list)
    return R_cmd

# Compare the concordance of two sets of labels 
def compute_df_concordance(y_ref, y_query, sort_list=True):
    """Compare the concordance of two sets of labels 

    Args:
        y_ref (list/array): reference label.
        y_query (list/array): query label.

    Returns:
        df_concordance (df): concordance matrix
    """
    list_ref = list(set(y_ref))
    list_query = list(set(y_query))
    if sort_list:
        list_ref.sort()
        list_query.sort()
    df_concordance = pd.DataFrame(index = list_query, columns = list_ref, data=0) 
    for val_query,val_ref in product(list_query, list_ref):
        df_concordance.loc[val_query, val_ref] = np.sum((y_ref==val_ref) &
                                                        (y_query==val_query))
    return df_concordance

# parse clustering annotation result
def parse_cluster_result(input_data, ref_list=[]):
    """Compare the concordance of two sets of labels 

    Args:
        input_data (adata, with clustering result): input clustered data 
        ref_list (list): list of annotations to match the cluster result with

    Returns:
        df_cluster_annotation (df): parsed clustering result
    """
    temp_adata = input_data.copy()
    temp_adata.obs['age_yo'] = ['young' if x in ['1m', '3m'] else 'old'
                                for x in temp_adata.obs['age']]
    # A sorted age list
    age_list = [int(x[:-1])for x in set(temp_adata.obs['age'])]
    age_list.sort()
    age_list = ['%dm'%x for x in age_list]
    # Cluster list 
    cluster_list = [str(x) for x in
                    np.sort(np.array(list(set(temp_adata.obs['leiden'])), dtype=int))]
    # Build cluster annotation
    df_cluster_annotation = pd.DataFrame(index=cluster_list)
    df_cluster_annotation['cluster_size'] = [np.sum(temp_adata.obs['leiden']==x)
                                             for x in cluster_list]
    # Add count for each age
    temp_df = compute_df_concordance(temp_adata.obs['leiden'], temp_adata.obs['age'])
    temp_df = temp_df.loc[age_list]
    temp_df.index = ['%s.ct'%x for x in age_list]
    df_cluster_annotation = df_cluster_annotation.join(temp_df.transpose())
    # Add normalized proportion for each age
    temp_df = temp_df.divide(temp_df.sum(axis=1), axis='rows')*1000
    temp_df = temp_df.divide(temp_df.sum(axis=0), axis='columns')
    temp_df.index = ['%s.prop'%x for x in age_list]
    df_cluster_annotation = df_cluster_annotation.join(temp_df.transpose())
    # Do the same for age_yo
    age_yo_list = ['young', 'old']
    temp_df = compute_df_concordance(temp_adata.obs['leiden'], temp_adata.obs['age_yo'])
    temp_df = temp_df.loc[age_yo_list]
    temp_df.index = ['%s.ct'%x for x in age_yo_list]
    df_cluster_annotation = df_cluster_annotation.join(temp_df.transpose())
    # Add normalized proportion for each age_yo
    temp_df = temp_df.divide(temp_df.sum(axis=1), axis='rows')*1000
    temp_df = temp_df.divide(temp_df.sum(axis=0), axis='columns')
    temp_df.index = ['%s.prop'%x for x in age_yo_list]
    df_cluster_annotation = df_cluster_annotation.join(temp_df.transpose())
    # Add ref annotation
    for ref in ref_list:
        if ref not in temp_adata.obs.keys():
            print('%s not valid'%ref)
            continue
        temp_df = compute_df_concordance(temp_adata.obs[ref], temp_adata.obs['leiden'])
        if 'nan' in temp_df.columns:
            df_cluster_annotation['%s.nan.prop'%ref] = temp_df['nan'] / temp_df.sum(axis=1)
            temp_df = temp_df.drop('nan', axis=1)
        else:
            df_cluster_annotation['%s.nan.prop'%ref] = 0
        temp_df = temp_df.divide(temp_df.sum(axis=1), axis='rows')
        if temp_df.shape[1]>0:
            df_cluster_annotation['%s.top.name'%ref] = \
                [temp_df.columns[np.argsort(temp_df.loc[x])[-1]] for x in cluster_list]
            df_cluster_annotation['%s.top.prop'%ref] = \
                [temp_df.loc[x, df_cluster_annotation['%s.top.name'%ref].loc[x]] 
                 for x in cluster_list]
        # print(cluster_list)
        if temp_df.shape[1]>1:
            df_cluster_annotation['%s.second.name'%ref] = \
                [temp_df.columns[np.argsort(temp_df.loc[x])[-2]] 
                 for x in cluster_list]
            # print(temp_df)
            df_cluster_annotation['%s.second.prop'%ref] = \
                [temp_df.loc[x, df_cluster_annotation['%s.second.name'%ref].loc[x]]
                 for x in cluster_list]
        if temp_df.shape[1]>2:
            df_cluster_annotation['%s.third.name'%ref] = \
                [temp_df.columns[np.argsort(temp_df.loc[x])[-3]] 
                 for x in cluster_list]
            # print(temp_df)
            df_cluster_annotation['%s.third.prop'%ref] = \
                [temp_df.loc[x, df_cluster_annotation['%s.third.name'%ref].loc[x]]
                 for x in cluster_list]
        temp = []
        for x in cluster_list:
            if (df_cluster_annotation.loc[x, '%s.nan.prop'%ref]>0.8) | \
                (df_cluster_annotation.loc[x, 'cluster_size']<20):
                temp.append('nan')
            else:
                if (df_cluster_annotation.loc[x, '%s.top.prop'%ref]>0.9):
                    temp.append('%s %d%%'%
                                (df_cluster_annotation['%s.top.name'%ref].loc[x],
                                 100*df_cluster_annotation['%s.top.prop'%ref].loc[x]))
                elif ((df_cluster_annotation.loc[x, ['%s.top.prop'%ref, '%s.second.prop'%ref]].sum()>0.9 )):
                    temp.append('%s %d%%/\n%s %d%%'%
                                (df_cluster_annotation['%s.top.name'%ref].loc[x],
                                 100*df_cluster_annotation['%s.top.prop'%ref].loc[x],
                                 df_cluster_annotation['%s.second.name'%ref].loc[x],
                                 100*df_cluster_annotation['%s.second.prop'%ref].loc[x]))
                elif ((df_cluster_annotation.loc[x, ['%s.top.prop'%ref, '%s.second.prop'%ref, '%s.third.prop'%ref]].sum()>0.9 )):
                    temp.append('%s %d%%/\n%s %d%%/\n%s %d%%'%
                                (df_cluster_annotation['%s.top.name'%ref].loc[x],
                                 100*df_cluster_annotation['%s.top.prop'%ref].loc[x],
                                 df_cluster_annotation['%s.second.name'%ref].loc[x],
                                 100*df_cluster_annotation['%s.second.prop'%ref].loc[x],
                                 df_cluster_annotation['%s.third.name'%ref].loc[x],
                                 100*df_cluster_annotation['%s.third.prop'%ref].loc[x]))
                else:
                    temp.append('nan')
        df_cluster_annotation['%s.name'%ref] = temp
    return df_cluster_annotation