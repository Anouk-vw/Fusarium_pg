"""
Read dupfile and create heatmap
"""
from sys import argv
import pandas as pd
from natsort import natsort_keygen, index_natsorted
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

#get gene location
def gp_to_dict(gp):
    gp_dict = {}
    with open(gp) as gp_file:
        for line in gp_file:
            line = line.strip()
            gene = line.split('\t')[0]
            chrom = line.split('\t')[1]
            start = line.split('\t')[3]
            stop = line.split('\t')[4]
            gp_dict[gene] = f'{chrom}-{start}-{stop}'
    return gp_dict
    
def MCgff_to_dict(MCgff):
    MCgff_dict = {}
    with open(MCgff) as gff_file:
        for line in gff_file:
            line = line.strip()
            gene = line.split('\t')[1]
            chrom = line.split('\t')[0]
            start = line.split('\t')[2]
            stop = line.split('\t')[3]
            MCgff_dict[gene] = f'{chrom}-{start}-{stop}'
    return MCgff_dict

#get gene type df
def get_dup_type(in_file):
    type_dict = {}
    with open(in_file) as dup_file:
        for line in dup_file:
            gene, gene_type = line.strip().split('\t')
            #if gene_type != '0':
            type_dict[gene] = gene_type
    return type_dict

def link_dict(dup_dict, loc_dict):
    linked_dict = {}
    for k,v in dup_dict.items():
        linked_dict[k] = loc_dict[k].split('-')
        linked_dict[k].extend(dup_dict[k])
    return linked_dict

def sort_df(dict_dups):
    df = pd.DataFrame(dict_dups)
    df = df.T
    df.columns = ['chrom', 'start','stop','type']
    df = df.sort_values(["chrom", "start"], key=natsort_keygen())
    return df
    
#create column per type
def col_pertype(df,col):
    for dtype in set(list(df[col])):
        df[f'col_{dtype}'] = (df[col] == dtype)
    return df

#cummulative position

#heatmap row per type

if __name__ == '__main__':
    gp_locdict = MCgff_to_dict(argv[1]) #gp_to_dict(argv[1])
    dup_gene_dict = get_dup_type(argv[2])
    ldict = link_dict(dup_gene_dict, gp_locdict)
    df_toplot = sort_df(ldict)
    col_type_df = col_pertype(df_toplot, 'type')
    #ax, fig = plt.subplots(figsize=(10,5))
    
    #adaptive_col = ((col_type_df['chrom'] == 'TR4_II5_Chr1') & \
    #        (col_type_df['stop'].astype(int) < 1800000)) | \
    #        (col_type_df['chrom'] == 'TR4_II5_Chr12')
    #row_col = ['orange' if x else 'black' for x in adaptive_col]
            
    print(col_type_df)
    cols = sorted([x for x in col_type_df.columns if 'col_' in x])
    print(cols)
    sns.clustermap(\
    col_type_df.loc[:,cols].T,\
    row_cluster=False, col_cluster=False,\
    figsize=(10,2),\
    cmap=ListedColormap(['white', '#2BB5C2']), xticklabels=False,\
    yticklabels = ['Singleton','Dispersed','Proximal','Tandem','Segmental'])#,\
    #col_colors=row_col)
    
    plt.tight_layout()
    plt.savefig(f'{argv[1].split("/")[-1]}.pdf')
    #print(col_type_df[col_type_df['col_4']])
