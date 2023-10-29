#!/usr/bin/env python
# coding: utf-8

# In[1]:

"""
script to go from orthofinder to gene age classification 
file: gene-name LCA-phyla age-number

python Gene_age.py ../Orthofinder_results/Results_Feb23_1/Orthogroups/Orthogroups.GeneCount.tsv ../Orthofinder_results/Results_Feb23_1/Orthogroups/Orthogroups.tsv 



"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from natsort import natsorted
from natsort import natsort_keygen
from sys import argv


def read_orthogroups(OrthoFinder_table):
    #'Results_Feb23_1/Orthogroups/Orthogroups.GeneCount.tsv'
    OG_counts = {}
    i = 0
    with open(OrthoFinder_table) as OG_counts_file:
        for line in OG_counts_file:
            #if i < 100:
                i+=1
                line=  line.strip()
                if line.split('\t')[0] == 'Orthogroup':
                    OG_counts[line.split('\t')[0]] = line.split('\t')[1:]
                else:
                    OG_counts[line.split('\t')[0]] = [int(x) for x in line.split('\t')[1:]]
    print(f"Number of orthogroups: {len(OG_counts)}\nNumber of genomes:{len(OG_counts['Orthogroup'])}")
    return OG_counts

def list_of_Phyla(Phyla_ordered, accessions_df, Fob_names):
    """Return list of IDs in OrthoFinder ordered on Phyla"""
    #get list of IDs ordered on Phlya (.loc orders the df)
    ordered_df = accessions_df.set_index('Phylum').loc[Phyla_ordered]
    
    #initialize col ''Phylum_toreturn''
    ordered_df['Phylum_toreturn'] = ordered_df.index
    #hypocreales (not phyla but ordered) for resolution
    ordered_df.loc[ordered_df['Order'] == 'hypocreales','Phylum_toreturn'] = 'hypocreales'
    print(ordered_df.loc[ordered_df['Order'] == 'hypocreales'])
    #FusOx (not phyla, needed for resolution)
    ordered_df.loc[ordered_df['Portal'] == 'Fusoxrap1.prot','Phylum_toreturn'] = 'FusOx'

    #get Phyla names
    ordered_phyl_col = list(ordered_df['Phylum_toreturn'])
    ordered_phyl_col.extend(['FusOx' for x in Fob_names])
    print(np.unique(ordered_phyl_col))
   
    #get genome IDs ordered on Phyla
    ordered_IDs = list(ordered_df['Portal'])

    ordered_PhylumIDs = [f'{x}.prot' for x in ordered_IDs]
    ordered_PhylumIDs.extend(Fob_names)
    return ordered_PhylumIDs, ordered_phyl_col

def plot_OrthoDF_wrows(ordered_Phyla):
    col_dict = {'blastocladiomycota':'black', 'chytridiomycota':'grey', \
    'zoopagomycota':'blue','mucoromycota':'green', 'basidiomycota':'yellow', 'ascomycota':'red'}
    row_colors = [col_dict[x] for x in ordered_Phyla]
    row_colors.extend(['orange' for x in Fob_names])
    sns.clustermap(OG_df.loc[ordered_Phylum, OGs_P_Fob], vmax = 1, row_colors=row_colors,row_cluster=False, col_cluster=False,)

def get_OGs_inX(OG_df, IDs_on_Phylum, x, phyla_col):
    """Get list of OGs with presence in genomes in X"""
    OGs_in_x = OG_df.loc[:,list((OG_df.loc[x] != 0).sum() != 0)].columns
     
    x_genes_df = OG_df.loc[IDs_on_Phylum, OGs_in_x]
    
    x_genes_df_binary = (x_genes_df != 0)
    x_genes_df_binary['Phylum'] = phyla_col
    
    return x_genes_df_binary

def fill_dict(Phyla_oldtonew, binary_df):
    #Change Phylum to All_fungi
    age_dict = {'blastocladiomycota':'all_fungi', 'chytridiomycota':'all_fungi', 'zoopagomycota':'all_fungi',\
                         'mucoromycota':'all_fungi', 'basidiomycota':'all_fungi', 'ascomycota':'ascomycota'}
    binary_df['Phylum'] = [age_dict[x] if x in age_dict else x for x in binary_df['Phylum']]
    OG_per_Phylum = binary_df.groupby('Phylum').sum()
    Phyl_dict = {}
    Phyl_dict_percentage = {}
    #per file in list (ordered old to new)
    for Phyl in Phyla_oldtonew:
        #get OGs in Phylum (!=0)
        OGs_Phyl = list(OG_per_Phylum.loc[:,OG_per_Phylum.loc[Phyl] != 0].columns)
        #Add Phylum to dict {Phlyum:OGs in Phylum}
        Phyl_dict[Phyl] = OGs_Phyl
        
        #total organisms from Phylum (changes of a OG to occur there)
        total_org = (binary_df['Phylum'] == Phyl).sum()
        #OG occuring in x of the organisms:
        OG_percentage = OG_per_Phylum.loc[Phyl] / total_org >= 0.2 #more than 20% Presence
        OG_percentage_list = list(OG_per_Phylum.columns[OG_percentage])
        
        Phyl_dict_percentage[Phyl] = OG_percentage_list
        print(f"{Phyl}\ttotal: {len(OG_per_Phylum.loc[:,OG_per_Phylum.loc[Phyl] != 0].columns)}\t     percentage:\t{len(OG_percentage_list)}")
    return Phyl_dict, Phyl_dict_percentage

def get_gene_names(Orthofinder_out):
    #'Results_Feb23_1/Orthogroups/Orthogroups.tsv'
    OG_gene_names = {}
    with open(Orthofinder_out) as OG_counts_file:
        for line_num, line in enumerate(OG_counts_file):
            if line_num == 0:
                header = line.split('\t')
            #if line_num > 0 and line_num < 10:
            else:
                OG_gene_names[line.split('\t')[0]] =  {header[i+1]:x for i,x in enumerate(line.split('\t')[1:])}
    return OG_gene_names
    
def get_genes_Phylum(OGs_phylum, OG_gene_dict, strain):
    genes = []
    for OG in OGs_phylum:
        genes.extend([i.replace(' ','') for i in OG_gene_dict[OG][strain].split(',')])
    return genes

def age_per_gene(Phyla_old_to_new, Phyla_dict, OG_gene_names):
    """
    Phyla_dict = {Phyla:OGs}
    """
    #Phyla_old_to_new.reverse()
    Present_in_older = []
    dict_LCA = {}
    dict_gene_phyl = {}
    for Phyl in Phyla_old_to_new:
        #get genes in this Phylum
        genes_in_phyl = get_genes_Phylum(Phyla_dict[Phyl], OG_gene_names, 'TR4_II5.new_protein')
        #remove genes present in older Phylum, to make sure the last common presence is retained
        genes_most_recent = set(genes_in_phyl) - set(Present_in_older)
        #append all genes occuring in this phyla, to make sure the gene occurs in oldest Phyla
        Present_in_older.extend(genes_in_phyl)
        
        for gene in genes_most_recent:
            dict_gene_phyl[gene] = Phyl
        
        dict_LCA[Phyl] = genes_most_recent
    return dict_LCA
    
#80 are lost due to percentage (low percentage in all?)
#others are presumably unique. TODO
#Fob_genes_df.sum(axis = 1)['TR4_II5.new_protein']

def wrap(Ortho_res, Phyla_ordered, accessions, Ortho_res_genes, Phyla_ordered_complete):
    OG_count_dict = read_orthogroups(Ortho_res)
    OG_OF_df = pd.DataFrame(OG_count_dict)
    OG_OF_df =OG_OF_df.set_index('Orthogroup')
    
    Fob_names = [x for x in OG_OF_df.index if x.endswith('.new_protein')]
    
    #IDs on Phyla
    Phyl_id, Phyl_col = list_of_Phyla(Phyla_ordered, accessions, Fob_names)
    
    OG_Fob_df_binary = get_OGs_inX(OG_OF_df, Phyl_id, Fob_names, Phyl_col)
    
    OGs_in_Phylum, OGs_in_perc_Phylum = fill_dict(Phyla_ordered_complete, OG_Fob_df_binary)
    
    #Phyla_ordered_complete = ['all_fungi', 'ascomycota', 'hypocreales','FusOx']
    map_gene2OG = get_gene_names(Ortho_res_genes)
    
    gene_per_Phyla_dict = age_per_gene(Phyla_ordered_complete, OGs_in_perc_Phylum, map_gene2OG)
    
    return gene_per_Phyla_dict

if __name__ == '__main__':
    accessions_df = pd.read_csv('../trace_genes/gene_age/genome_accession_published_subset.csv')
    Phyla_old_to_new = ['blastocladiomycota', 'chytridiomycota', 'zoopagomycota','mucoromycota', 'basidiomycota', 'ascomycota']
    #Phyla_old_to_new_complete = ['blastocladiomycota', 'chytridiomycota', 'zoopagomycota','mucoromycota', 'basidiomycota', 'ascomycota', 'hypocreales','FusOx']
    Phyla_old_to_new_complete = ['all_fungi', 'ascomycota', 'hypocreales','FusOx']
    gene_age_d = wrap(argv[1], Phyla_old_to_new, accessions_df, argv[2], Phyla_old_to_new_complete)
    
    with open('gene_age_20perc_allfungi.txt', 'w+') as output:
        for k,v in gene_age_d.items():
            for gene in v:
                if len(gene) > 1:
                    output.write(f'{gene}\t{k}\t{Phyla_old_to_new_complete.index(k)}\n')
