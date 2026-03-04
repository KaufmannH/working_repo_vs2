import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import openpyxl
import seaborn as sns
import matplotlib.pyplot as plt



# select only 3 months old mice
def select_3m():
    print("Starting to load raw anndata object.")
    adata = sc.read_h5ad('/home/hkaufm49/working_repo/TMS/pansci/data/GSE247719_PanSci_Myeloid_cell_adata.h5ad')
    print("loaded adata object.")
    pansci_3m_wt = adata[ (adata.obs['Age_group']== "03_months") & (adata.obs['Genotype']== "WT")]
    pansci_3m_wt.write("/home/hkaufm49/working_repo/TMS/pansci/data/pansci_3m_wt.h5ad")
    print("Saved only 3 months old WT mice ")
    return(pansci_3m_wt)

#df_3m_wt = select_3m()

def save_cell_type_list():
    adata = sc.read_h5ad("/home/hkaufm49/working_repo/TMS/pansci/data/pansci_3m_wt.h5ad")
    sub_cell_type_list = pd.DataFrame(adata.obs['Sub_cell_type'].unique(), columns=['Sub_cell_type'])
    sub_cell_type_list.to_excel("/home/hkaufm49/working_repo/TMS/pansci/data/unique_sub_cell_types.xlsx", index=False)
    print("Saved cell type list.")

#save_cell_type_list()


def prep_pansci_df():
    adata = sc.read_h5ad("/home/hkaufm49/working_repo/TMS/pansci/data/pansci_3m_wt.h5ad")
    adata.obs['sex'] = adata.obs['Sex'].str.lower()
    adata.obs['tissue'] = adata.obs['Sub_cell_type'].str.split('-').str[-1].str.lower()
    adata.obs['age'] = adata.obs['Age_group'].str.split('_').str[0].astype(int)
    adata.obs['cell_type'] = adata.obs['Immune_subtype'].str.lower().str.replace(" ", "_")
    # generate cluster numbers
    unique_combos = adata.obs[['tissue', 'cell_type']].drop_duplicates().reset_index(drop=True)
    cluster_numbers = []
    for tissue in unique_combos['tissue'].unique():
        tissue_combos = unique_combos[unique_combos['tissue'] == tissue]
        cluster_numbers.extend(range(1, len(tissue_combos) + 1))
    unique_combos['cluster_number'] = cluster_numbers
    adata.obs = adata.obs.merge(unique_combos, on=['tissue', 'cell_type'], how='left')
    adata.obs = adata.obs.sort_values(by=['tissue', 'cell_type', 'cluster_number'])
    #adata.obs[['tissue', 'cell_type', 'cluster_number']].drop_duplicates().head(40)
    adata.obs['cluster_id'] =  adata.obs['tissue'].astype(str) + "_" +  adata.obs['age'].astype(str) + "m_" + adata.obs['sex'].astype(str) + "_" + adata.obs['cluster_number'].astype(str)
    adata.obs['mouse_id'] =  adata.obs['age'].astype(str) + "_" + adata.obs['ID'].astype(str) + "_" + adata.obs['Sex'].str[0]
    adata.obs = adata.obs[["cluster_id" , "cell_type", "tissue", "age", "mouse_id", "Gene_count"]]
    adata.write("/home/hkaufm49/working_repo/TMS/pansci/data/pansci_preped.h5ad")
    return(adata)

#preped_df = prep_pansci_df()


def select_cells():
    adata = sc.read_h5ad("/home/hkaufm49/working_repo/TMS/pansci/data/pansci_preped.h5ad")

    narrowed_cell_type_list = ['chil3+_alveolar_macrophages', 
                            'inflammatory_monocytes',
                            'patrolling_monocytes',
                            'col14a1+_Itga8+_macrophages',
                            'proliferating_macrophages',
                            'lyve1+_colec12+_macrophages',
                            'type-2_conventional_dendritic cells',
                            'type-1_conventional_dendritic cells',
                            'migratory_dendritic_cells',
                            'retnla+_cd226+_macrophages',
                            'cfd+_scd1+_macrophages',
                            'intestinal_macrophages',
                            'gpc6+_hspg2+_macrophages',
                            'kupffer_cells',
                            'mmp12+_mmp19+_macrophages',
                            'muscle-resident_macrophages',
                            'colq+_gpm6b+_macrophages',
                            'fcgr4+_itgad+_macrophages']
    adata2 = adata[adata.obs['cell_type'].isin(narrowed_cell_type_list)].copy()
    adata2.write("/home/hkaufm49/working_repo/TMS/pansci/data/pansci_innate_cells.h5ad")
    return(adata2)
#selected_cells = select_cells()


def downsample():
    adata = sc.read_h5ad("/home/hkaufm49/working_repo/TMS/pansci/data/pansci_innate_cells.h5ad")
    n_cells_to_sample = 5000
    np.random.seed(42)
    random_indices = np.random.choice(adata.obs.index, size=n_cells_to_sample, replace=False)
    adata_small = adata[random_indices].copy()
    print(adata_small.obs)
    adata_small.write("/home/hkaufm49/working_repo/TMS/pansci/data/pansci_downsampled.h5ad")
    return(adata_small)
#downsampled = downsample()



def gene_per_cell_cutoff(adata):
    adata = sc.read_h5ad("/home/hkaufm49/working_repo/TMS/pansci/data/pansci_innate_cells.h5ad")
    adata_small = adata[adata.obs['Gene_count'] > 500]
    adata_small.write("/home/hkaufm49/working_repo/TMS/pansci/data/pansci_500cutoff.h5ad")
    return(adata_small)

#small_adata = gene_per_cell_cutoff(selected_cells)





# not working: a lot of means are 0? 
def calc_mean(): 
    pansci_cleaned = sc.read_h5ad("/home/hkaufm49/working_repo/TMS/pansci/data/pansci_500cutoff.h5ad")
    
    X = pansci_cleaned.X
    cluster_ids = pansci_cleaned.obs['cluster_id'].to_numpy()
    gene_names = pansci_cleaned.var['gene_name']


    from scipy.sparse import issparse
    if issparse(X):
        X = X.tocoo()

    long_df = pd.DataFrame({
        'cluster_id': cluster_ids[X.row],
        'gene': gene_names[X.col],
        'expression': X.data})
    long_df.head()
    plt.figure(figsize=(10, 6))
    sns.histplot(data=long_df, x='expression', kde=True, bins=10)
    plt.title('Distribution of Gene Expression')
    plt.xlabel('Expression Level')
    plt.ylabel('Frequency')
    plt.savefig("/home/hkaufm49/working_repo/TMS/pansci/plots/hist_test.png")


    mean_expression_df = (
        long_df
        .groupby(['cluster_id', 'gene'], as_index=True)
        .agg(gmean=('expression', 'mean')))
    mean_expression_df.reset_index()

    # add metadata
    metadata_df = pansci_cleaned.obs[['cluster_id',  'cell_type', 'tissue', 'age', 'mouse_id', 'Gene_count']].drop_duplicates()
    # Extract metadata from AnnData
    # Extract metadata from AnnData
    mean_expression_df = mean_expression_df.merge(metadata_df, on='cluster_id', how='left')


    mean_expression_df.head()
    mean_expression_df['gmean'].describe()
    mean_expression_df['gmean'].max()

    # reporting
    plt.figure(figsize=(10, 6))
    sns.histplot(long_df['expression'])
    plt.xlabel('Mean Gene Expression')
    plt.ylabel('Number of Gene-Cluster Combinations')
    plt.title('Distribution of Mean Gene Expression across Clusters')
    plt.savefig("/home/hkaufm49/working_repo/TMS/pansci/plots/hist_test.png")







# not working: a lot of means are 0? 
def calc_gmean_old():
    pansci_cleaned = sc.read_h5ad("/home/hkaufm49/working_repo/TMS/pansci/data/pansci_500cutoff.h5ad")


    #pansci_cleaned = small_adata

    X = pansci_cleaned.X

    # avoid dense format
    X = X.tocoo()  # COO format (row, col, data)
    cell_ids = pansci_cleaned.obs_names.to_numpy()
    gene_names = pansci_cleaned.var['gene_name'].to_numpy()

    long_df = pd.DataFrame({
        'cell_id': cell_ids[X.row],
        'gene': gene_names[X.col],
        'expression': X.data
    })
    long_df.head()
    #long_df_t = (
       # long_df_filterd
        #.groupby(['cell_id'], as_index=False)
       # .agg({'gene': 'count'}))

    # make df with other metadata info
    obs_df = pansci_cleaned.obs.reset_index().rename(columns={'sample': 'cell_id'})
    
    merged_df = long_df.merge(obs_df, on=['cell_id'], how='left')
    #merged_df.to_csv("/home/hkaufm49/working_repo/TMS/pansci/data/pansci_df_long_merged.csv")


    mean_expression_df = (
        merged_df
        .groupby(['gene', 'cluster_id'], as_index=False)
        .agg({
            'expression': 'mean',
            'cell_type' : 'first',
            'age': 'first',
            'tissue': 'first'})
        .rename(columns={'expression': 'gmean'}))

    mean_expression_df.head()
    mean_expression_df.to_csv("/home/hkaufm49/working_repo/TMS/pansci/data/pansci_df_gmean.csv")

   
    mean_expression_df['gmean'].describe()
    mean_expression_df['gmean'].max()

    # reporting
    sns.histplot(mean_expression_df, binwidth=1)
    plt.xlabel('Mean expression (gmean)')
    plt.ylabel('Number of genes')
    plt.savefig("/home/hkaufm49/working_repo/TMS/pansci/plots/hist_test2.png")


    pandf = pd.read_csv("/home/hkaufm49/working_repo/TMS/pansci/data/pansci_df_gmean.csv")
    print(pandf)

