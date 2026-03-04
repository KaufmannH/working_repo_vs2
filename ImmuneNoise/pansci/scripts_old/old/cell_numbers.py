import seaborn as sns
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


adata = sc.read_h5ad("/home/hkaufm49/working_repo/TMS/pansci/data/pansci_preped.h5ad")

def cell_type_dist():

    df_cell_types = adata.obs[['cell_type']].copy()
    perc_cell_type = (
        df_cell_types['cell_type']
        .value_counts()
        .reset_index())

    perc_cell_type['percent'] = perc_cell_type['count'] / perc_cell_type['count'].sum() * 100
    perc_cell_type = perc_cell_type.sort_values('percent', ascending=False).reset_index(drop=True)
    order = perc_cell_type['cell_type'].tolist()


    plt.figure(figsize=(12, 7))

    ax = sns.barplot(
        data=perc_cell_type,
        x='cell_type', y='percent',
        order=order, 
        palette= sns.color_palette("husl", len(perc_cell_type)))
    
    for bar, count in zip(ax.patches, perc_cell_type['count']):
        height = bar.get_height()
        colour = bar.get_facecolor()

        ax.text(
            bar.get_x() + bar.get_width() / 2,
            height + 1,
            str(int(count)),
            ha='center', va='bottom', fontsize=9,
            rotation=45,
            color=colour)

    ax.set_ylim(0, perc_cell_type['percent'].max() + 10)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_ylabel("Percentage of all cells")
    ax.set_xlabel("Cell type")
    ax.set_title("Cell type composition")
    ax.yaxis.set_major_formatter(plt.matplotlib.ticker.PercentFormatter(100))
    sns.despine()
    plt.tight_layout()
    plt.tight_layout()
    plt.savefig("/home/hkaufm49/working_repo/TMS/pansci/plots/cell_type_composition.png", dpi=300)
    plt.close()

#cell_type_dist()





def organ_cell_dist():

    df_organs = adata.obs[['tissue']].copy()

    perc_organs = (
        df_organs['tissue']
        .value_counts()
        .reset_index())
    perc_organs['percent'] = perc_organs['count'] / perc_organs['count'].sum() * 100
    perc_organs = perc_organs.sort_values('percent', ascending=False).reset_index(drop=True)
    order = perc_organs['tissue'].tolist()


    plt.figure(figsize=(12, 7))

    ax = sns.barplot(
        data=perc_organs,
        x='tissue', y='percent',
        order=order, 
        palette= sns.color_palette("husl", len(perc_organs)))

    for bar, count in zip(ax.patches, perc_organs['count']):
        height = bar.get_height()
        colour = bar.get_facecolor()

        ax.text(
            bar.get_x() + bar.get_width() / 2,
            height + 1,
            str(int(count)),
            ha='center', va='bottom', fontsize=9,
            rotation=45,
            color=colour)

    ax.set_ylim(0, perc_organs['percent'].max() + 10)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_ylabel("Percentage of all cells")
    ax.set_xlabel("Tissue")
    ax.set_title("Tissue cell composition")
    ax.yaxis.set_major_formatter(plt.matplotlib.ticker.PercentFormatter(100))
    sns.despine()
    plt.tight_layout()
    plt.savefig("/home/hkaufm49/working_repo/TMS/pansci/plots/organ_composition.png", dpi=300)
    plt.close()

organ_cell_dist()