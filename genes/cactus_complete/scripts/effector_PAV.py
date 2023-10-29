"""
script to parse brocoli (+effP+busco)
Save tables
"""

from sys import argv
import pandas as pd
import random
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import sys

def get_OG_cat(dataframe, cutoff):
    """
    For dataframe, 
    """
    cat =[]
    counts = []
    genes_c = []
    for i,x in enumerate(list((dataframe!= 0).sum(axis = 1))):
        x = int(x)
        counts.append(x)
        #sum each row (total number of genes) and list. 
        #From this list get element i (row travesed by loop)
        genes_c.append(int(list(dataframe.sum(axis = 1))[i]))
        if int(x) >= cutoff:
            cat.append('core')
        if x < cutoff and x >= cutoff-10:
            cat.append('softcore')
        if x < cutoff-10 and x>=2:
            cat.append('accessory')
        if x < 2:
            cat.append('unique')
    dataframe['category'] = cat
    dataframe['counts'] = counts
    #representing number of genes in all genomes
    dataframe['gene_counts'] = genes_c
    return dataframe

def genes_OG(b_table_names):
    #get all genes per OG
    genes_OG = {}
    with open(b_table_names) as broc_names:
        for line in broc_names:
            if line.startswith('#'):
                continue
            else:
                OG = line.split('\t')[0]
                ##This returns a list if one genome has more 
                ##than one gene in the cluster. Flatten list
                genes_OG[OG] = list(np.concatenate([gene.split(' ') \
                for gene in line.strip().split('\t')[1:] if len(gene) > 1]))
    return genes_OG

def effector_list(effector_names):
    """
    list all effectors in eff_file
    """
    eff_list = []
    with open(effector_names) as eff_file:
        for line in eff_file:
            eff_list.append(line.strip())
    return set(eff_list)

def is_gene_effector(gene_name, effector_list):
    return gene_name in effector_list

def OG_effectors(OG_genes_dict, effector_list):
    """List all OGs with effectors"""
    eff_dict = {}
    effs = []
    for OG,genes in OG_genes_dict.items():
        for gene in genes:
            eff = is_gene_effector(str(gene), effector_list)
            try:
                eff_dict[OG].append(eff)
            except KeyError:
                eff_dict[OG] = [eff]
    return eff_dict

def get_unique(unique, effector_list):
    """
    Get unique eff from unique file
    """
    eff_in_u = []
    u_eff = {}
    with open(unique) as u:
        for line in u:
            if line.strip() in effector_list:
                u_eff[line.strip()] = 1
            else:
                u_eff[line.strip()] = 0
    return u_eff

def get_percentage(df, column, column_dev):
    percentage_effectors = df[column].astype(int) / df[column_dev].astype(int)
    percentage_effectors[percentage_effectors != 0.0]
    print(f'no {column}:', len(percentage_effectors[percentage_effectors == 0.0]),
      f'less than 50% {column}:', len(percentage_effectors[(percentage_effectors < 0.5) & (percentage_effectors > 0.0)]),
      f'more than 50% {column}:', len(percentage_effectors[percentage_effectors >= 0.5]))
    return percentage_effectors[percentage_effectors >= 0.5].index

def fill_gene_df(counts, names, eff, uniques,busco):
    """
    Take all dataframes and counts, append to dataframe
    """
    #get counts
    tableCounts = pd.read_csv(counts, sep = '\t', index_col=0)
    #add category and count (based on tableCounts)
    catdf = get_OG_cat(tableCounts, 69)
    #create Orthogroup dictonary {OG:[genes]}
    OG_dict = genes_OG(names)
    #Read effector list
    print(eff)
    effectors = effector_list(eff)
    busco_genes = effector_list(busco)
    #make dictionary {OG:[boolean_list, effector_true/false]}
    OG_eff_dict = OG_effectors(OG_dict, effectors)
    OG_busco_dict = OG_effectors(OG_dict, busco_genes)
    #make dictionary {OG:total_effectors}
    OG_eff_dict = {k:sum(v) for k,v in OG_eff_dict.items()}
    OG_busco_dict = {k:sum(v) for k,v in OG_busco_dict.items()}
    #read uniques and get {unique:effector_True/False}
    u_dict = get_unique(uniques, effectors)
    u_dict_busco = get_unique(uniques, busco_genes)
    #make dataframe with uniques
    udf = pd.DataFrame.from_dict(u_dict, orient = 'index')
    udf_busco = pd.DataFrame.from_dict(u_dict_busco, orient = 'index')
    #make dataframe Orthogroups [OG][eff_count]
    odf = pd.DataFrame.from_dict(OG_eff_dict, orient = 'index')
    busco_df = pd.DataFrame.from_dict(OG_busco_dict, orient = 'index')
    #concatenate unique dict and orthogroup dict (add unique's as index)
    concat = pd.concat([udf, odf])
    conc_busco = pd.concat([udf_busco, busco_df])
    #Combine category, count and effector_count in dataframe [OG-cat-count-eff_count]
    comb_df = pd.concat([tableCounts[['category', 'counts', 'gene_counts']], concat, conc_busco], axis = 1)
    #replace Nan of unique values (appended later, not in counts df)
    comb_df['category'] = comb_df['category'].replace(np.nan, 'Unique')
    comb_df['counts'] = comb_df['counts'].replace(np.nan, 1)
    comb_df['gene_counts'] = comb_df['gene_counts'].replace(np.nan, 1)
    comb_df.columns = ['category', 'counts', 'gene_counts', 'effectors', 'busco']
    eff_50_OGs = get_percentage(comb_df, 'effectors', 'gene_counts')
    comb_df['eff_perc'] = ['effector' if x else 'no_effector' \
        for x in comb_df.index.isin(eff_50_OGs)]
    busco_50_OGs = get_percentage(comb_df, 'busco', 'gene_counts')
    comb_df['busco_perc'] = ['busco' if x else 'no_busco' \
        for x in comb_df.index.isin(busco_50_OGs)]
    return comb_df

def get_unique_perc(genome, singleton_dir):
    """
    get unique genes and determine effector percentage
    """
    u_genes = []
    with open(f"{singleton_dir}/{genome}_singletons") as s:
        for line in s:
            u_genes.append(line.strip())
    uniques = [x for x in u_genes if x in effectors]
    perc_unique_eff = (len(uniques)/len(u_genes))
    return perc_unique_eff
    
def per_genome_eff_class(genome, read_df, effectors, extra_df, singleton_dir, ext):
    """
    return percentage effectors per category for genome
    """
    OG_dataframe = read_df[read_df[f"{genome}{ext}"].notna()]
    TR4_OGs = list(OG_dataframe.index)
    if genome == '15_KTG06-B2':
        effectors_TR4 = [x for x in effectors if '15_KTG06-B' in x]
    else:
        effectors_TR4 = [x.strip() for x in effectors if genome in x]
    TR4_genes_OG = OG_dataframe[OG_dataframe[f'{genome}{ext}'].isin(effectors_TR4)].index
    not_unique = []
    [not_unique.extend(y) for y in [x.split(' ') for x in OG_dataframe[f'{genome}{ext}']]]
    eff_OGs = list(TR4_genes_OG)
    uniques = get_unique_perc(genome, singleton_dir)
    perc_gene_is_eff = np.divide(np.unique(list(extra_df.loc[eff_OGs]['category']), return_counts = True)[1], np.unique(list(df.loc[TR4_OGs]['category']), return_counts = True)[1]) * 100
    if len(perc_gene_is_eff) == 3:
        perc_gene_is_eff = np.append(perc_gene_is_eff, perc_unique_eff)
    return perc_gene_is_eff

def eff_content(names_file, counts_file, effector_file, singleton_dir, ext):
    """
    Return list of effector content per location per genome
    """
    df = fill_gene_df(counts_file, names, eff)
    read_df = pd.read_csv(names_file, sep = '\t', index_col = 0)
    with open(eff) as eff_file:
        effector_list_ = [x.rstrip() for x in eff_file.readlines()]
    for g in [x.replace(ext, '') for x in read_df.columns][1:-1]:
        genomes_order.append(g)
        acc, core, sc, u = per_genome_eff_class(g, read_df, effector_list_, df, singleton_dir, ext)
        accessory_effector.append(acc)
        core_effector.append(core)
        sc_effector.append(sc)
        u_effector.append(u)
    return core_effector, sc_effector, accessory_effector, u_effector

def plot_dist_effector(names_file, counts_file, effector_file):
    """
    plot distribution of percentage of effectors
    """
    eff_contents = eff_content(names_file, counts_file, effector_file)
    ax = sns.boxplot(data=eff_contents, color = "0.75")
    ax.set_ylabel("Percentage Effectors", labelpad=10)
    plt.xticks([0,1,2,3], ['core', 'softcore', 'accessory', 'unique'])
    sns.stripplot(data=[core_effector, sc_effector, accessory_effector, u_effector],  color="0", alpha = 0.5)
    #plt.savefig("boxplot_effectordistribution.pdf")

def get_clean_table(counts_file, cutoff, ext):
    table = get_OG_cat(pd.read_csv(counts_file, sep = '\t', index_col=0), cutoff)
    table.columns = [x.split(ext)[0] for x in table.columns]
    to_plot = table.iloc[:,0:69]
    to_plot = to_plot > 1
    table = to_plot
    return table

def plot_PAV(counts_file, names_file, eff_file, cutoff, ordering, out):
    """
    Plot PAV of effectors either all or only accessory
    """
    table_cluster = get_clean_table(counts_file)
    table_cluster.iloc[:,0:69]
    df = fill_gene_df(counts_file, names_file, eff_file)
    accessory_effectors = df[(df[0] == 1) & (df['category'] != 'core') & (df['category'] != 'Unique')]
    all_effectors = df[(df[0] == 1)]
    ordered = table_cluster.reindex(list(phyl_order), axis = 1)
    sns.clustermap(ordered.reindex(accessory_effectors.index), \
    figsize = (30,20), col_cluster = False, \
    cmap=ListedColormap(['white', '#B30050']), cbar_pos = None)
    plt.savefig(out)
    return

def get_specific_OG(counts_file, meta, meta_col, meta_val, cutoff):
    """
    does not work
    counts = argv[1]
    cf = 69
    meta = pd.read_csv('/home/anouk/SPP_FOCReSequencingMetadata.csv', sep = ',')
    Tr4_specific = get_specific_OG(counts, meta, 'Tr4.R1', 'TR4', cf)
    print(Tr4_specific)
    """
    df = get_clean_table(counts_file, cutoff)
    True_in_names = [x for x in list(meta[meta[meta_col] == meta_val]['IsolateCode'])]
    False_in_names = [x for x in list(meta[meta[meta_col] != meta_val]['IsolateCode'])]
    index_true = df[list(df.reindex(columns = True_in_names).all(axis = 1))]
    index_true = index_true[index_true.reindex(columns = False_in_names).sum(axis = 1) < 1]
    return index_true.index

def genes_from_OG(genome, names_df, gene_list):
    TR4_genes_OG = names_df[names_df[genome].notna()][genome]
    acc_genes_TR4 = TR4_genes_OG.reindex(gene_list)[TR4_genes_OG.reindex(gene_list).notna()]
    return acc_genes_TR4

if __name__ == '__main__':
    #counts = argv[1]
    cf = 69
    meta = pd.read_csv('/home/anouk/SPP_FOCReSequencingMetadata.csv', sep = ',')
    #Tr4_specific = get_specific_OG(counts, meta, 'Tr4.R1', 'TR4', cf)
    #print(Tr4_specific)
    genes_path  = \
    '/home/anouk/anouk2/pangenome/cactus_complete/genes/'
    counts = \
    f'{genes_path}/brocoli_results/dir_step3/table_OGs_protein_counts.txt'
    names = \
    f'{genes_path}/brocoli_results/dir_step3/table_OGs_protein_names.txt'
    eff = \
    f'{genes_path}/Predicted_effectors_noMINI.txt'
    uniques =\
    f'{genes_path}/brocoli_results/unique_genes.txt'
    busco_genes = \
    f'{genes_path}/new_protein_files/busco_annotation/all_complete_busco_genes.txt'
    #read broccoli
    #table = get_OG_cat(pd.read_csv(counts, sep = '\t', index_col=0), 69)
    #table.columns = [x.replace(f'{ext}', '') for x in table.columns]
    info_df = fill_gene_df(counts, names, eff, uniques, busco_genes)
    info_df.to_csv(argv[1], header = ['category', 'species_count', \
            'gene_count','effector','busco', 'eff_perc', 'busco_perc'])
