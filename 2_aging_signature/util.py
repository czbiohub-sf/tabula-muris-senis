# import scanpy.api as sc
import scanpy as sc
import pandas as pd
import numpy as np
import scipy as sp
# import scipy.stats as stats
import time
import statsmodels.formula.api as smf
import statsmodels.api as sm
import os
from anndata import read_h5ad
from anndata import AnnData
from itertools import product
import matplotlib.pyplot as plt
from sys import getsizeof

def label_bar(v_p, v_x, v_y, method='fwer', alpha=0.05, dh=.05):
    
    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh *= (ax_y1 - ax_y0)
    
    if method == 'fwer':
        v_p = v_p * v_p.shape[0]
    v_y = v_y.clip(min=0)
        
    for i in np.arange(v_p.shape[0]):
        text='*' if v_p[i]<alpha else 'n.s.'
        plt.annotate(text, [v_x[i], v_y[i]+dh], ha='center', va='bottom')
        
    plt.ylim([ax_y0, max(ax_y1, v_y.max()+dh+(ax_y1 - ax_y0)*0.2)])

def label_bardif(text, lx, rx, ly, ry, 
                 dh=.05, barh=.05, loc='upper',
                 fs=None, maxasterix=None):
    """ 
    Annotate barplot with p-values.

    :param num1: number of left bar to put bracket over
    :param num2: number of right bar to put bracket over
    :param data: string to write or number for generating asterixes
    :param center: centers of all bars (like plt.bar() input)
    :param height: heights of all bars (like plt.bar() input)
    :param yerr: yerrs of all bars (like plt.bar() input)
    :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
    :param barh: bar height in axes coordinates (0 to 1)
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    """
    
    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh *= (ax_y1 - ax_y0)
    barh *= (ax_y1 - ax_y0)
    
    if loc=='upper':
        y = max(ly, ry) + dh

        barx = [lx, lx, rx, rx]
        bary = [y, y+barh, y+barh, y]
        mid = ((lx+rx)/2, y+1.1*barh)

        plt.plot(barx, bary, c='black')
        plt.ylim([ax_y0, max(ax_y1, (y+barh)*1.2)])

        kwargs = dict(ha='center', va='bottom')
        if fs is not None:
            kwargs['fontsize'] = fs
            
        plt.text(*mid, text, **kwargs)
        
    else:
        y = max(ly, ry) - dh

        barx = [lx, lx, rx, rx]
        bary = [y, y-barh, y-barh, y]
        mid = ((lx+rx)/2, y-1.1*barh)

        plt.plot(barx, bary, c='black')
        plt.ylim([min(ax_y0, (y-barh)*1.2), ax_y1])

        kwargs = dict(ha='center', va='top')
        if fs is not None:
            kwargs['fontsize'] = fs
            
        plt.text(*mid, text, **kwargs)
        
    
def test_overlap(list1, list2, list_background):
    
    set1 = set(list1)
    set2 = set(list2)
    set_background = set(list_background)
    
    n1 = len(set1)
    n2 = len(set2)
    n_overlap = len(set1 & set2)
    n_other = len(set_background - set1 - set2)
    
    oddsratio, pvalue = sp.stats.fisher_exact([[n_other, n1-n_overlap],
                                               [n2-n_overlap, n_overlap]])
    return oddsratio,pvalue

def test_twosample(v1, v2):
    
    t,p = scipy.stats.ttest_ind(v1, v2, equal_var=True)
    return t,p

def get_p_two_point(mean1, se1, mean2=0, se2=0):
    dif_mean = (mean1-mean2)
    dif_se = np.sqrt(se1**2 + se2**2)
    z_score = dif_mean/dif_se
    p_val = 2*sp.stats.norm.cdf(-np.absolute(z_score))
    return p_val

def get_p_two_point_v(v_mean1, v_se1, v_mean2=0, v_se2=0):
    v_dif_mean = (v_mean1-v_mean2)
    v_dif_se = np.sqrt(v_se1**2 + v_se2**2)
    v_z_score = v_dif_mean/v_dif_se
    v_p_val = 2*sp.stats.norm.cdf(-np.absolute(v_z_score))
    return v_p_val

def sizeof(a):
    return getsizeof(a)/1024/1024

def nested_partition(df_data, response, regressor_list, covariate_list=None):
    
    # Regress out covariates 
    if covariate_list is not None:
        pass
    
    # 
    df_res = pd.DataFrame(index=['rsquared', 'rsquared_adj'], 
                          columns=regressor_list,
                          data=0)
    
    # Add regressors one-by-one
    v_y = df_data[response].values
    temp_list = []
    for term in regressor_list:
        temp_list.append(term)
        mat_X = df_data[temp_list].values
        mat_X = sm.add_constant(mat_X, prepend=False)
        res = sm.OLS(v_y, mat_X).fit()
        df_res.loc['rsquared', term] = res.rsquared
        df_res.loc['rsquared_adj', term] = res.rsquared_adj
    
    return df_res
    

def get_sharing(input_df_coef, input_df_fdr, 
                coef_thres=0.005, fdr_thres=0.01, verbose=False):
    
    df_coef = input_df_coef.copy()
    df_fdr = input_df_fdr.copy()
    
    n_gene,n_cond = df_coef.shape
    cond_list = list(df_coef.columns)
    df_sharing = pd.DataFrame(index=cond_list, columns=cond_list, data=0)
    
    if verbose: 
        print('# get_sharing: n_cond=%d, n_gene=%d'
              %(n_cond, n_gene))
    
    for i_cond1 in np.arange(n_cond):
        
        cond1 = cond_list[i_cond1]
        v_coef1 = df_coef[cond1].values
        v_fdr1 = df_fdr[cond1].values
        v_sig1 = (v_fdr1<fdr_thres) & (np.absolute(v_coef1)>coef_thres)
        
        for i_cond2 in np.arange(i_cond1, n_cond):
            
            if i_cond1==i_cond2: 
                prop_sharing = 0 
                
            else:
                cond2 = cond_list[i_cond2]
                v_coef2 = df_coef[cond2].values
                v_fdr2 = df_fdr[cond2].values
                v_sig2 = (v_fdr2<fdr_thres) & (np.absolute(v_coef2)>coef_thres)
                
                v_coef2[v_coef2==0] = 1e12
                r_coef = v_coef1/v_coef2
                v_coef_overlap = (r_coef>0.5) & (r_coef<2)
                
                n_sig = (v_sig1 | v_sig2).sum()
                
                n_sharing = (v_sig1 & v_sig2 & v_coef_overlap).sum()
                
                prop_sharing = n_sharing/n_sig
                                
            df_sharing.iloc[i_cond1,i_cond2] = prop_sharing
            df_sharing.iloc[i_cond2,i_cond1] = prop_sharing
                
    return df_sharing

def compute_aging_score_model(gene_list_up, gene_list_down):
    
    gene_list =  gene_list_up + gene_list_down
    df_model = pd.DataFrame(index = gene_list, columns=['coef'], data=0)
    
    df_model.loc[gene_list_up, 'coef'] = 1/ len(gene_list_up)
    df_model.loc[gene_list_down, 'coef'] = -1/len(gene_list_down)
        
    return df_model

def compute_aging_score_model_qw(adata, gene_list_up, gene_list_down):
    
    gene_list =  gene_list_up + gene_list_down
    df_model = pd.DataFrame(index = gene_list, columns=['coef'], data=0)
    
    temp_X = adata[:, gene_list].X.toarray()
    v_range = np.quantile(temp_X, 0.95, axis=0) - np.quantile(temp_X, 0.05, axis=0)
    df_model['coef'] = 1/ v_range.clip(min=0.1)
    df_model.loc[gene_list_up, 'coef'] = df_model.loc[gene_list_up, 'coef'] \
                                            / df_model.loc[gene_list_up, 'coef'].sum()
    df_model.loc[gene_list_down, 'coef'] = -df_model.loc[gene_list_down, 'coef'] \
                                            / df_model.loc[gene_list_down, 'coef'].sum()

    return df_model 

# def compute_aging_score_model(adata, gene_list_up, gene_list_down, option='equal_weight', verbose=False):
#     
#     if option not in ['equal_weight', 'quantile_weight']: 
#         print('# compute_aging_score_model: option %s not defined, exit with None'%option)
#         return None
#     
#     gene_list =  gene_list_up + gene_list_down
#     df_model = pd.DataFrame(index = gene_list, columns=['coef'], data=0)
#     
#     if option=='equal_weight':
#         df_model.loc[gene_list_up, 'coef'] = 1/ len(gene_list_up)
#         df_model.loc[gene_list_down, 'coef'] = -1/ len(gene_list_down)
#         
#     if option=='quantile_weight':
#         temp_X = adata[:, gene_list].X.todense()
#         v_range = np.quantile(temp_X, 0.95, axis=0) - np.quantile(temp_X, 0.05, axis=0)
#         df_model['coef'] = 1/ v_range.clip(min=0.5)
#         df_model.loc[gene_list_up, 'coef'] = df_model.loc[gene_list_up, 'coef'] \
#                                                 / np.sum(df_model.loc[gene_list_up, 'coef'])
#         df_model.loc[gene_list_down, 'coef'] = -df_model.loc[gene_list_down, 'coef'] \
#                                                 / np.sum(df_model.loc[gene_list_down, 'coef'])
#         
#     return df_model

def compute_aging_score(adata,
                        df_model,
                        flag_correct_background=True,
                        verbose=False):
        
    if verbose: 
        start_time=time.time()
    
    # Set parameters
    n_gene_model = df_model.shape[0]
    gene_list = list(set(adata.var_names) & set(df_model.index))
    df_model = df_model.loc[gene_list].copy()
    
    # Rescale the weights 
    ind_select = df_model['coef']>0
    total_coef_up = np.absolute(df_model.loc[ind_select, 'coef'].sum())
    df_model.loc[ind_select, 'coef'] = df_model.loc[ind_select, 'coef'] / total_coef_up
    ind_select = df_model['coef']<0
    total_coef_down = np.absolute(df_model.loc[ind_select, 'coef'].sum())
    df_model.loc[ind_select, 'coef'] = df_model.loc[ind_select, 'coef'] / total_coef_down
    
    v_weight = df_model['coef'].values
    
    if verbose:
        print('# compute_aging_score: n_cell=%d, n_gene_adata=%d, n_gene_model=%d, n_gene_overlap=%d'%
              (adata.shape[0], adata.shape[1], n_gene_model, len(gene_list)))
        print('# compute_aging_score: per_cell_overall_exp_for_model_genes=%0.1f'%
              (adata[:, gene_list].X.sum(axis=1).mean()))
        print('# compute_aging_score: total_coef_up=%0.2f, total_coef_down=%0.2f'%
              (total_coef_up, total_coef_down))
    
    # Compute aging scores 
    df_aging_score = pd.DataFrame(index=adata.obs.index, columns=['score'], data=0)
    temp_X = adata[:, gene_list].X
    v_weight = v_weight.reshape([-1, 1])
    df_aging_score['score'] = (temp_X.dot(v_weight))[:,0]
    
    if flag_correct_background:
        v_mean,v_var = get_sparse_var(adata.X, axis=1)
        v_std = np.sqrt(v_var)
        df_aging_score['score'] = (df_aging_score['score'] - v_mean*v_weight.sum()) / \
                                    (v_std * np.sqrt((v_weight**2).sum()))
        
    # Covariates
    df_aging_score = df_aging_score.join(adata.obs[['age', 'age_num', 'sex', 'analyte']])
    
    if verbose: 
        print('# compute_aging_score: finished, time=%0.1fs'%(time.time()-start_time))
    
    return df_aging_score

def get_sparse_var(sparse_X, axis=0):
    
    if sp.sparse.issparse(sparse_X):
        v_mean = sparse_X.mean(axis=axis)
        v_mean = np.array(v_mean).reshape([-1])
        v_var = sparse_X.power(2).mean(axis=axis)
        v_var = np.array(v_var).reshape([-1])
        v_var = v_var - v_mean**2
    else:
        v_mean = sparse_X.mean(axis=axis)
        v_var = sparse_X.var(axis=axis)
    
    return v_mean,v_var

# def compute_aging_score(adata, df_model, verbose=False):
#     
#     # fixit: correct for background 
#     
#     if verbose: 
#         start_time=time.time()
#     
#     # Set parameters
#     n_gene_model = df_model.shape[0]
#     gene_list = list(set(adata.var_names) & set(df_model.index))
#     df_model = df_model.loc[gene_list].copy()
#     
#     # Rescale the weights 
#     ind_select = df_model['coef']>0
#     total_coef_up = np.absolute(df_model.loc[ind_select, 'coef'].sum())
#     df_model.loc[ind_select, 'coef'] = df_model.loc[ind_select, 'coef'] / total_coef_up
#     ind_select = df_model['coef']<0
#     total_coef_down = np.absolute(df_model.loc[ind_select, 'coef'].sum())
#     df_model.loc[ind_select, 'coef'] = df_model.loc[ind_select, 'coef'] / total_coef_down
#     
#     v_weight = df_model['coef'].values
#     
#     if verbose:
#         print('# compute_aging_score: n_cell=%d, n_gene_adata=%d, n_gene_model=%d, n_gene_overlap=%d'%
#               (adata.shape[0], adata.shape[1], n_gene_model, len(gene_list)))
#         print('# compute_aging_score: per_cell_overall_exp_for_model_genes=%0.1f'%
#               (adata[:, gene_list].X.sum(axis=1).mean()))
#         print('# compute_aging_score: total_coef_up=%0.2f, total_coef_down=%0.2f'%
#               (total_coef_up, total_coef_down))
#     
#     # Compute aging scores 
#     df_aging_score = pd.DataFrame(index=adata.obs.index, columns=['score'], data=0)
#     # df_aging_score = pd.DataFrame(index=adata.obs.index, 
#     #                               columns=['score', 'score.sex_adj', 'score.age_sex_adj'],
#     #                               data=0)
#     temp_X = adata[:, gene_list].X
#     v_weight = v_weight.reshape([-1, 1])
#     df_aging_score['score'] = (temp_X.dot(v_weight))[:,0]
#     
#     # Compute adjusted score
#     df_aging_score = df_aging_score.join(adata.obs[['age', 'age_num', 'sex', 'analyte']])
#     
#     # # Compute sex adjusted score
#     # result = smf.ols(formula="score ~ sex", data=df_aging_score).fit()
#     # df_aging_score['score.sex_adj'] = result.resid
#     # 
#     # # Compute age-and-sex adjusted score
#     # result = smf.ols(formula="score ~ age + sex + age*sex", data=df_aging_score).fit()
#     # # result = sm.ols(formula="score ~ age_num + sex + age_num*sex", data=df_aging_score).fit()
#     # df_aging_score['score.age_sex_adj'] = result.resid
#     
#     if verbose: 
#         print('# compute_aging_score: finished, time=%0.1fs'%(time.time()-start_time))
#     
#     return df_aging_score

def meta_analysis(effects, se, method='random', weights=None):
    # From Omer Weissbrod
    assert method in ['fixed', 'random']
    d = effects
    variances = se**2
    
    #compute random-effects variance tau2
    vwts = 1.0 / variances
    fixedsumm = vwts.dot(d) / vwts.sum()    
    Q = np.sum(((d - fixedsumm)**2) / variances)
    df = len(d)-1
    tau2 = np.maximum(0, (Q-df) / (vwts.sum() - vwts.dot(vwts) / vwts.sum()))
    
    #defing weights
    if weights is None:
        if method == 'fixed':
            wt = 1.0 / variances
        else:
            wt = 1.0 / (variances + tau2)
    else:
        wt = weights
    
    #compute summtest
    summ = wt.dot(d) / wt.sum()
    if method == 'fixed':
        varsum = np.sum(wt*wt*variances) / (np.sum(wt)**2)
    else:
        varsum = np.sum(wt*wt*(variances+tau2)) / (np.sum(wt)**2)
    ###summtest = summ / np.sqrt(varsum)
    
    summary=summ
    se_summary=np.sqrt(varsum)
    
    return summary, se_summary

def load_DGE_res(file_path, dname='bulk.tissue', version='1e4'):
    
    """Load DGE results 

    Args:
        file_path (str): file path. Should contain both FACS data facs_filtered.h5ad and droplet data droplet_filtered.h5ad
        
        version: one of '1e4' or 'tpm'
            '1e4': DGE results correponding to the normalization that total_ct_per_cell=1e4
            'tpm': DGE results correponding to the normalization that total_ct_per_cell=1e6

    Returns:
        adata_combine (AnnData): Combined data for FACS and droplet
    """
        
    if dname=='bulk.tissue':
        RESULT_PATH = file_path + '/DGE_result/DGE_bulk.%s'%version
    elif dname=='facs.tissue':
        RESULT_PATH = file_path + '/DGE_result/DGE_facs_tissue.%s'%version
    elif dname=='facs_old.tissue':
        RESULT_PATH = file_path + '/DGE_result/DGE_facs_old_tissue'
    elif dname=='facs.tc':
        # RESULT_PATH = file_path + '/DGE_result/DGE_facs_tc'
        RESULT_PATH = file_path + '/DGE_result/DGE_facs_tc.%s'%version
    elif dname=='facs_old.tc':
        RESULT_PATH = file_path + '/DGE_result/DGE_facs_old_tc'
    elif dname=='droplet.tissue':
        RESULT_PATH = file_path + '/DGE_result/DGE_droplet_tissue.%s'%version
    elif dname=='droplet_old.tissue':
        RESULT_PATH = file_path + '/DGE_result/DGE_droplet_old_tissue'
    elif dname=='droplet.tc':
        # RESULT_PATH = file_path + '/DGE_result/DGE_droplet_tc'
        RESULT_PATH = file_path + '/DGE_result/DGE_droplet_tc.%s'%version
    elif dname=='droplet_old.tc':
        RESULT_PATH = file_path + '/DGE_result/DGE_droplet_old_tc'
    else:
        return None,None
    
    method_name,tc_name = dname.split('.')
    df_info = pd.read_csv(file_path + '/DGE_result/%s.%s_info'%(method_name,tc_name), sep=' ', index_col=0)
        
    dic_res = {}
    dir_list = os.listdir(RESULT_PATH)
    for i_analyte,analyte in enumerate(df_info.index):
        analyte = analyte.replace(' ', '_')
        fname = '%s.csv'%analyte
        if fname not in dir_list:
            print('%d, %s missing'%(i_analyte+1,fname))
        else:
            dic_res[analyte] = pd.read_csv(RESULT_PATH+'/'+fname, sep=',', index_col=0)
            
    # fixit: change analyte name
    
    return df_info,dic_res

# # Code to get the raw_data_combined.h5ad for ma_data
# import os
# import scipy.io
# import anndata
# 
# data_path = "/n/groups/price/martin/tms_gene_data/Ma_Cell_2020_data"
# file_list = os.listdir(data_path)
# file_list = ['%s_%s'%(x.split('_')[0], x.split('_')[1]) for x in file_list
#              if 'mtx' in x]
# 
# dic_tissue = {'Aorta':'Aorta', 'BM':'Marrow', 'BAT':'BAT', 'Brain':'Brain', 'Liver':'Liver',
#               'Muscle':'Limb_Muscle', 'Skin':'Skin', 'WAT':'WAT', 'Kidney':'Kidney'}
# dic_sex = {'F':'female', 'M':'male'}
# dic_age = {'Y':'3m', 'O':'27m', 'CR':'27m'}
# dic_cond = {'Y':'AL', 'O':'AL', 'CR':'CR'}
# 
# adata = None
# for fname in file_list:
#     
#     prefix = data_path + '/%s'%fname
#     # Parse fname
#     tissue = dic_tissue[fname.split('_')[1].split('-')[0]]
#     sex = dic_sex[fname.split('_')[1].split('-')[1]]
#     age = dic_age[fname.split('_')[1].split('-')[2]]
#     cond = dic_cond[fname.split('_')[1].split('-')[2]]
#     print(fname, tissue, sex, age, cond)
#     
#     mat = scipy.io.mmread(prefix+'_matrix.mtx.gz')
#     mat = sp.sparse.csr_matrix(mat)
#     barcodes_path = prefix+"_barcodes.tsv.gz"
#     
#     df_gene = pd.read_csv(prefix+"_genes.tsv.gz", sep='\t', compression='gzip', header=None)
#     df_gene.columns = ['GENE_ID', 'GENE']
#     df_gene.index = df_gene['GENE']
#     
#     df_obs = pd.read_csv(prefix+"_barcodes.tsv.gz", sep='\t', compression='gzip', header=None)
#     df_obs.columns = ['BARCODE']
#     df_obs.index = df_obs['BARCODE']
#     
#     df_obs['tissue'] = tissue
#     df_obs['sex'] = sex
#     df_obs['age'] = age
#     df_obs['cond'] = cond
#     
#     temp_adata = anndata.AnnData(mat.T, obs=df_obs, var=df_gene)
#     temp_adata.var_names_make_unique()
#     
#     if adata is None:
#         adata = temp_adata.copy()
#     else:
#         adata = adata.concatenate(temp_adata)
#         
# adata.write(data_path+'/raw_data_combined.h5ad')

# # Code to get the filtered_data.h5ad for ma_data
# data_path = "/n/groups/price/martin/tms_gene_data/Ma_Cell_2020_data"
# adata = read_h5ad(data_path+'/raw_data_combined.h5ad')
# 
# sc.pp.filter_genes(adata, min_cells=5)
# sc.pp.filter_cells(adata, min_genes=500)
# 
# adata.obs['n_counts'] = adata.X.sum(axis=1)
# adata.obs['n_genes'] = (adata.X>0).sum(axis=1)
# adata.obs['age_num'] = [int(x.replace('m','')) for x in adata.obs['age']]
# mt_gene_list = [x for x in adata.var_names if 'Mt-' in x]
# print(mt_gene_list)
# adata.obs['mt_counts'] = adata[:,mt_gene_list].X.sum(axis=1)
# adata.obs['mt_frac'] = adata.obs['mt_counts']/adata.obs['n_counts']
# adata = adata[adata.obs['mt_frac']<0.1,:]
# 
# adata.write(data_path+'/filtered_data.h5ad')

def load_ma_data(root_data_path, 
                     flag_size_factor=True, total_ct_per_cell=1e4, 
                     flag_log1p=True):
    
    """Load normalized data from Kimmel et al, GR, 2019
    1. Size factor normalization to counts per 1 million (total_ct_per_cell)
    2. log(x+1) transform

    Args:
        file_path (str): file path. Should contain ./Kimmel_GR_2019_data

    Returns:
        adata_combine (AnnData): Combined data for kidney, lung, and spleen
    """
    
    # Load filtered data
    file_path=root_data_path+'/Ma_Cell_2020_data' 
    adata = read_h5ad(file_path + '/filtered_data.h5ad')
    
    # Size factor normalization
    if flag_size_factor == True:
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=total_ct_per_cell)
        
    # log(x+1) transform
    if flag_log1p == True:
        sc.pp.log1p(adata)
    return adata



# Code used to curated the data for load_kowalczyk_data
# root_data_path=DATA_PATH
# file_path=root_data_path+'/Kowalczyk_GR_2015' 
# # df_data = pd.read_excel(file_path+'/GSE59114_C57BL6_GEO_all.xlsx')
# df_data = pd.read_csv(file_path+'/GSE59114_C57BL6_GEO_all.txt',
#                       sep='\t', skiprows=1, index_col=0, header=0)
# cell_list = [x for x in df_data if x.split('_')[0] in ['young', 'old']]
# df_data = df_data[cell_list].copy().T
# df_data.columns = [x.replace("'", '') for x in df_data.columns]
# 
# def get_celltype(x):
#     if 'LT_HSC' in x:
#         return 'LT_HSC'
#     if 'ST_HSC' in x:
#         return 'ST_HSC'
#     if 'MPP' in x:
#         return 'MPP'
#     return 'unknown'
# 
# df_obs = pd.DataFrame(index=df_data.index)
# df_obs['age'] = ['3m' if 'young' in x else '22m' for x in df_obs.index]
# df_obs['cell_type'] = [get_celltype(x) for x in df_obs.index]
# 
# import anndata
# adata = anndata.AnnData(df_data, obs=df_obs)
# adata.X = scipy.sparse.csr_matrix(adata.X)
# adata.write(file_path+'/GSE59114_C57BL6_GEO_all.h5ad')

def load_kowalczyk_data(root_data_path, total_ct_per_cell=1e4):
    
    """Load normalized data from Kowalczyk et al, GR, 2015
    The data is normalized with log2(tpm+1) 
    Fixit: change to 1e4

    Args:
        file_path (str): file path. Should contain ./Kowalczyk_GR_2015/GSE59114_C57BL6_GEO_all.h5ad

    Returns:
        adata (AnnData): 
    """
    
    # read file 
    file_path=root_data_path+'/Kowalczyk_GR_2015' 
    adata=read_h5ad(file_path+'/GSE59114_C57BL6_GEO_all.h5ad')
    adata.obs['age_num'] = [int(x.replace('m','')) for x in adata.obs['age']]
    adata.obs['sex'] = 'female'
    adata.var_names_make_unique(join='_')
    
    if total_ct_per_cell!=1e6:
        temp_X = adata.X.todense()
        temp_X = np.exp2(temp_X) - 1
        temp_X = temp_X*(total_ct_per_cell/1e6)
        adata.X = sp.sparse.csr_matrix(temp_X)
        sc.pp.log1p(adata)
    
    return adata

def load_kimmel_data(root_data_path, 
                     flag_size_factor=True, total_ct_per_cell=1e4, 
                     flag_log1p=True):
    
    """Load normalized data from Kimmel et al, GR, 2019
    1. Size factor normalization to counts per 1 million (total_ct_per_cell)
    2. log(x+1) transform

    Args:
        file_path (str): file path. Should contain ./Kimmel_GR_2019_data

    Returns:
        adata_combine (AnnData): Combined data for kidney, lung, and spleen
    """
    
    # Load filtered data
    file_path=root_data_path+'/Kimmel_GR_2019_data' 
    adata_kidney = read_h5ad(file_path + '/kidney.h5ad')
    adata_lung = read_h5ad(file_path + '/lung.h5ad')
    adata_spleen = read_h5ad(file_path + '/spleen.h5ad')
    
    # Size factor normalization
    if flag_size_factor == True:
        sc.pp.normalize_per_cell(adata_kidney, counts_per_cell_after=total_ct_per_cell)
        sc.pp.normalize_per_cell(adata_lung, counts_per_cell_after=total_ct_per_cell)
        sc.pp.normalize_per_cell(adata_spleen, counts_per_cell_after=total_ct_per_cell)
        
    # log(x+1) transform
    if flag_log1p == True:
        sc.pp.log1p(adata_kidney)
        sc.pp.log1p(adata_lung)
        sc.pp.log1p(adata_spleen)
    
    # Combine data 
    adata = adata_kidney.concatenate(adata_lung, adata_spleen,
                                     batch_key='batch_combine', join='inner')
    adata.obs['tissue'] = ''
    adata.obs.loc[adata.obs['batch_combine']=='0', 'tissue'] = 'Kidney'
    adata.obs.loc[adata.obs['batch_combine']=='1', 'tissue'] = 'Lung'
    adata.obs.loc[adata.obs['batch_combine']=='2', 'tissue'] = 'Spleen'
    
    adata.obs['sex'] = 'male'
    adata.obs['age_old'] = adata.obs['age'].values.copy()
    adata.obs['age'] = ['7m' if x=='young' else '22m' for x in adata.obs['age_old']]
    adata.obs['age_num'] = [7 if x=='young' else 22 for x in adata.obs['age_old']]
    
    return adata


# def load_normalized_data_by_tissue(file_path, data_name='facs.Aorta',
#                                    flag_size_factor=True, flag_log1p=True):
#     
#     data_name,tissue = data_name.split('.')
#     file_name = f"{data_name}.raw.{tissue}.h5ad"
#         
#     # Load filtered data
#     adata = read_h5ad(f'{file_path}/raw_adata_by_tissue/{file_name}')
#     
#     # Size factor normalization
#     if flag_size_factor == True:
#         sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
#         
#     # log(x+1) transform
#     if flag_log1p == True:
#         sc.pp.log1p(adata)
#         
#     return adata

def load_normalized_data(file_path, data_name='facs', 
                         flag_size_factor=True, 
                         total_ct_per_cell=1e4, 
                         flag_log1p=True):
    
    """load normalized data
    1. Load filtered data for both FACS and droplet
    2. Size factor normalization to counts per 1 million (total_ct_per_cell)
    3. log(x+1) transform
    4. Combine the data 

    Args:
        file_path (str): file path. Should contain both FACS data facs_filtered.h5ad and droplet data droplet_filtered.h5ad

    Returns:
        adata_combine (AnnData): Combined data for FACS and droplet
    """
        
    if data_name=='facs':
        file_name = 'tabula-muris-senis-facs-official-raw-obj.h5ad'
    elif data_name=='droplet':
        file_name = 'tabula-muris-senis-droplet-official-raw-obj.h5ad'
    elif data_name=='facs_old':
        file_name = 'facs_filtered.h5ad'
    elif data_name=='droplet_old':
        file_name = 'droplet_filtered.h5ad'
    else:
        return None
        
    # Load filtered data
    adata = read_h5ad(f'{file_path}/{file_name}')
    
    # Update annotations
    adata.obs['n_genes'] = (adata.X>0).sum(axis=1)
    adata.obs['n_counts'] = (adata.X).sum(axis=1)
    adata.obs['age_num'] = [int(x.replace('m','')) for x in adata.obs['age']]
    
    # Size factor normalization
    if flag_size_factor == True:
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=total_ct_per_cell)
        
    # log(x+1) transform
    if flag_log1p == True:
        sc.pp.log1p(adata)
        
    # Filter data
    if 'facs' in data_name:
        ind_select = adata.obs['age'].isin(['3m', '18m', '24m'])
        adata = adata[ind_select,]
    
    return adata

def load_normalized_data_bulk(file_path, flag_size_factor=True, total_ct_per_cell=1e4, flag_log1p=True):
    """load normalized bulk data
    1. Load filtered data for bulk RNA-Seq
    2. Size factor normalization to counts per 10 thousand
    3. log(x+1) transform

    Args:
        file_path (str): file path. Should contain the bulk data 190304_maca_bulk.h5ad
        log1p (bool): If to perform log1o transform

    Returns:
        adata_combine (AnnData): Combined data for FACS and droplet
    """
    # Load filtered data
    adata_bulk = read_h5ad(file_path+'/190304_maca_bulk.h5ad')
    adata_bulk.raw = adata_bulk
    # QC filtering
    sc.pp.filter_genes(adata_bulk, min_cells=3)
    sc.pp.filter_cells(adata_bulk, min_genes=250)
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
    if flag_size_factor:
        sc.pp.normalize_per_cell(adata_bulk, counts_per_cell_after=total_ct_per_cell)
    # log(x+1) transform
    if flag_log1p:
        sc.pp.log1p(adata_bulk)
    return adata_bulk