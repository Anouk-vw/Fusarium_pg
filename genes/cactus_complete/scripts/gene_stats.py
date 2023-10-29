"""
get busco list: 
grep 'Complete' *_prot_busco/run_hypocreales_odb10/full_table.tsv | cut -f 3 > all_complete_busco_genes.txt 
"""

import sys
sys.path.append('/home/anouk/anouk2/pangenome/cactus_complete/genes/')
import scripts.MIMPs_analysis as MIMP
import scripts.effector_PAV as eff_PAV
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from sys import argv
from Bio.Phylo.PAML import codeml
import scipy as sp
import itertools

# collect gene catergories: BUSCO, Core, Accessory, Effectors

def read_input(counts_table, names_table, cutoff, ext): 
	"""
	cutoff = Int. Occurences to consider core 
	ext = string to replace. Can be ''
	"""
	table = eff_PAV.get_OG_cat(pd.read_csv(counts_table, sep = '\t', index_col=0), cutoff)
	table.columns = [x.replace(f'{ext}', '') for x in table.columns]
	return table

def subset_categories(df, columns):
    group_OGs = {}
    for x in columns:
        for cat in set(list(df[x])):
            #0:100 to speed up plotoutline temporarily
            group_OGs[cat] = list(df[df[x] == cat].index)
    return group_OGs

def plot_distribution(OG_table, cutoff,binsize,save = False):
	"""
	Return histogram, distribution of gene PAV
	cutoff = Int. Max number of occurences
	binsize = Int. Size of bins
	save = Bool, TRUE to save plot 
	"""
	#range(2 to exclude singletons
	ax1 = plt.hist(OG_table['counts'], bins = range(2,cutoff,binsize))
	plt.rc('axes', axisbelow=True)
	plt.grid(axis = 'y', zorder = 0)
	plt.xlabel('Number of genomes')
	plt.ylabel('Number of OGs')
	if save:
		plt.savefig('hist_genePAV.pdf')
	else:
		plt.show()
	#print numbers:
	print(f"Core OGs:{len(OG_table[OG_table['category'] == 'core'])}\n\
	Softcore OGs:{len(OG_table[OG_table['category'] == 'softcore'])}\n\
	Accessory OGs:{len(OG_table[OG_table['category'] == 'accessory'])}\n\
	Unique OGs:{len(OG_table[OG_table['category'] == 'Unique'])}")
	return


#Number of effectors per category + MannwithenyU

#Length per category
def fasta_length(ID_OG, prefix, end):
    """
    Get length from fasta file for OG
    Return mean length of OG
    """
    lengths = []
    with open(f"{prefix}{ID_OG}{end}") as FastaFile:
        for name, seq in SimpleFastaParser(FastaFile):
            seqLen = len(seq)
            lengths.append(seqLen)
    return lengths

def toplot_length_per(grouped_OG_dict, prefix, end):
    """
    df = info df
    column = column to select on
    broccoli_df, 'category'
    """
    lists_toplot = []
    labels_plot = []
    for k,v in grouped_OG_dict.items():
        mean_len_OGs = [np.mean(fasta_length(OG, prefix, end)) for OG in v]
        lists_toplot.append(mean_len_OGs)
        labels_plot.append(k)
    
    return lists_toplot, labels_plot

def plot_box(list_of_list, labels_list, title, save = False):
    flat_list = []
    [flat_list.extend(x) for x in list_of_list]
    mean_all = np.median(flat_list)
    plt.subplots(figsize = (5,10))
    box = plt.boxplot(list_of_list, showfliers = False, \
            patch_artist=True, widths = 0.9)
    #box = plt.violinplot(list_of_list)
    colors = ['#146152', '#FFEC5C', '#44803F', '#B4CF66', '#FF5A33']
    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
    plt.xticks(range(1,len(list_of_list)+1), labels_list, \
            fontsize = 24,rotation = 90)
    plt.title(title)
    plt.ylabel(title, fontsize = 42)
    plt.hlines(mean_all, 0, len(list_of_list)+0.7)
    plt.yticks(fontsize = 24)
    plt.xlim(0.3,len(list_of_list)+0.7)
    if save:
        plt.savefig(f"boxplot_{title}.pdf", bbox_inches='tight')
        plt.close()
    else:
        plt.tight_layout()
        plt.show()
    """
    ax = sns.boxplot(data=list_of_list, color = "0.75")
    ax.set_ylabel("Count", labelpad=10)
    plt.xticks(range(0,len(list_of_list)), labels_list)
    sns.stripplot(data=list_of_list,  color="0", alpha = 0.5)
    plt.hlines(mean_all, 0, len(list_of_list)+1)
    plt.ylabel(title)
    if save:
        plt.savefig(f"boxplot_{title}_stripe.pdf")
        plt.close()
    else:
        plt.show()
    """
    return flat_list

#effector 'enrichment'
def toplot_eff_box(grouped_OGs, in_df):
    list_of_perc_eff = []
    labels_plot = []
    for k,v in grouped_OGs.items():
        #loop over list and get effectors per category
        num_effectors_OG = [in_df.loc[OG, 'effector'] for OG in v]
        num_genes_OG = [in_df.loc[OG, 'gene_count'] for OG in v]
        list_of_perc_eff.append(sum(num_effectors_OG) / sum(num_genes_OG))
        labels_plot.append(k)
    return list_of_perc_eff,labels_plot
    
#function -> perc annotated domains

#enrichment functions (eggnog)

#dn/ds
def parse_codeML(codeML_in):
    omega =[]
    try:
        results = codeml.read(codeML_in)
        keys =  results['pairwise'].keys()
        for kcomp in keys:
            subdict = results['pairwise'][kcomp]
            for k in keys:
                if k != kcomp:
                    #ds < 1, omega < 999
                    if float(subdict[k]['omega']) < 999:
                        omega.append(float(subdict[k]['omega']))
    except ValueError:
        #print(codeML_in)
        return False
    except KeyError:
        #print(codeML_in)
        return False
    return omega#.values()

def loop_codeml(in_dir):
    results_dict = {}
    for in_file in os.listdir(in_dir):
        omega_list = parse_codeML(f'{in_dir}{in_file}')
        if omega_list:
            OG = '_'.join(in_file.split('_')[0:2])
            results_dict[OG.replace('.txt','')] = np.median(omega_list)
    return results_dict

def toplot_dnds(OG_dndsdict, OG_group):
    list_of_dnds = []
    labels_plot = []
    for cat,OGs_list in OG_group.items():
        list_of_dnds.append([float(x) for O,x in OG_dndsdict.items()\
            if O in OGs_list])
        labels_plot.append(cat)
    return list_of_dnds,labels_plot

#TPM 
def parse_kallisto(input_file):
    RNA_dict = {}
    with open(input_file) as kallisto_file:
        for line in kallisto_file:
            gene, length, eff_length, count, tpm = line.strip().split('\t')
            RNA_dict[gene] = [length, eff_length, count, tpm]
    return RNA_dict

def map_OG_TR4(broc_dir):
    OG2TR4 = {}
    with open(f'{broc_dir}/dir_step3/table_OGs_protein_names.txt') as names:
        for line in names:
            for x in line.split('\t'):
                if 'OG' in x:
                    OG = x
                if 'TR4' in x:
                    TR4_gene = x
            OG2TR4[OG] = TR4_gene
            #for gene in TR4_gene.split(' '):
              #  OG2TR4[gene] = OG
    return OG2TR4
  
def get_mean_TPM_per_OG(kallisto_path, condition_list, OG_list, OG_TR4_dict):
    RNA_complete_tmp = []
    for cond in condition_list:
        print(f"{kallisto_path}{cond}")
        TPM_dict = {}
        for repl in ['A', 'B', 'C']:
            Kallisto_dict = parse_kallisto(f"{kallisto_path}{cond}.{repl}/abundance.tsv")
            for OG in OG_list:
                try:
                    tpm = Kallisto_dict[OG_TR4_dict[OG]][-1]
                except KeyError:
                    tpm = np.NaN
                try:
                    TPM_dict[OG].append(float(tpm))
                except KeyError:
                    TPM_dict[OG] = [float(tpm)]
        TPM_dict = {k:np.mean(v) for k,v in TPM_dict.items()}
        RNA_complete_tmp.append(TPM_dict)
    RNA_complete = {}
    for d in RNA_complete_tmp:
        for k,v in d.items():
            try:
                RNA_complete[k].append(v)
            except KeyError:
                RNA_complete[k] = [v]
    return RNA_complete
    
def toplot_tpm_per(grouped_OG_dict, TPM_table):
    """
    df = info df
    column = column to select on
    broccoli_df, 'category'
    """
    lists_toplot = []
    labels_plot = []
    for k,v in grouped_OG_dict.items():
        #0 = invitro, 1= 4dpi, 2=8dpi, 3 = 30dpi
        tpm_OGs = [TPM_table.loc[OG, 0] for OG in v]
        lists_toplot.append([x for x in tpm_OGs if x != 0])
        labels_plot.append(k)
    return lists_toplot, labels_plot
#prox to TE

#prox to other genes

#option to look only at subset of OGs

if __name__ == '__main__':
    brocoli_path  = '/home/anouk/anouk2/pangenome/cactus_complete/genes/brocoli_results'
    counts = f'{brocoli_path}/dir_step3/table_OGs_protein_counts.txt'
    names = f'{brocoli_path}/dir_step3/table_OGs_protein_names.txt'
    eff = '/home/anouk/anouk2/pangenome/cactus_complete/genes/Predicted_effectors_noMINI.txt'
    uniques = '/home/anouk/anouk2/pangenome/cactus_complete/genes/brocoli_results/unique_genes.txt'
    busco_genes = '/home/anouk/anouk2/pangenome/cactus_complete/genes/new_protein_files/busco_annotation/all_complete_busco_genes.txt'
    broc_table = read_input(counts, names, 69, '.new_protein.fasta')
    #For now ignore Unique
    info_df = pd.read_csv('../info_broccoli.csv', index_col = 0)
    info_df = info_df[info_df['category'] != 'Unique']

    #plot dist False = do not save
    ##plot_distribution(broc_table, 69, 1, False)
    
    #group_OGs:
    grouped_dict = \
        subset_categories(info_df, ['busco_perc', 'eff_perc', 'category'])
    grouped_dict.pop('no_effector')
    grouped_dict.pop('no_busco')
    print(grouped_dict.keys())

    #plot lengths
    print('plotting lengths')
    len_list, len_labels = toplot_length_per(grouped_dict,\
    '/home/anouk/anouk2/pangenome/cactus_complete/genes/align_nostop/',\
     '_prot.fasta_adapted.aln')
    list_reordered, label_list_reordered = [], []
    for k in ['busco', 'effector', 'core', 'softcore', 'accessory']:
        index = len_labels.index(k)
        list_reordered.append(len_list[index])
        label_list_reordered.append(k)
    flat_len = plot_box(list_reordered, label_list_reordered, 'length_violin', True)
    for i1,i2 in itertools.combinations(range(len(list_reordered)), 2):
        print('Length mannwithney U between',
                label_list_reordered[i1], f'mean: {np.mean(list_reordered[i1])}', \
                label_list_reordered[i2], f'mean: {np.mean(list_reordered[i2])}',
                sp.stats.mannwhitneyu(list_reordered[i1], list_reordered[i2]))
    [print(f"between {label_list_reordered[i]} {np.median(l)} and all {np.median(flat_len)}: \
    {sp.stats.mannwhitneyu(l, flat_len)}") for i,l in enumerate(list_reordered)]
    
    #print eff percentage
    print('calculating effectors')
    grouped_dicteff = \
        subset_categories(info_df, ['category', 'busco_perc'])
    grouped_dicteff.pop('no_busco')
    
    eff_list, eff_labels = toplot_eff_box(grouped_dicteff, info_df)
    print(eff_list, eff_labels)

    #plot TPM
    print('plotting TPM')
    OG2TR4_dict = map_OG_TR4(brocoli_path)
    TPM_dict = get_mean_TPM_per_OG\
    ('/home/anouk/anouk2/pangenome/cactus_complete/genes/kallisto/II5_kallisto/II5', \
            ['_Invitro', '_Infected4dpi', '_Infected8dpi', '_Infected30dpi'], broc_table.index, OG2TR4_dict)
    TPM_df = pd.DataFrame(TPM_dict).T.fillna(0)
    print(TPM_df)
    table_TPM = broc_table.join(TPM_df)
    list_tpm, labels_tpm = toplot_tpm_per(grouped_dict, table_TPM)
    list_reordered, label_list_reordered = [], []
    for k in ['busco', 'effector', 'core', 'softcore', 'accessory']:
        index = labels_tpm.index(k)
        list_reordered.append(list_tpm[index])
        label_list_reordered.append(k)
    flat_len = plot_box(list_reordered, label_list_reordered, 'TPM_Invitro_violin', True)
    for i1,i2 in itertools.combinations(range(len(list_reordered)), 2):
        print('Length mannwithney U between',
                label_list_reordered[i1], f'mean: {np.mean(list_reordered[i1])}', \
                label_list_reordered[i2], f'mean: {np.mean(list_reordered[i2])}',
                sp.stats.mannwhitneyu(list_reordered[i1], list_reordered[i2]))
    [print(f"between {label_list_reordered[i]} {np.median(l)} and all {np.median(flat_len)}: \
    {sp.stats.mannwhitneyu(l, flat_len)}") for i,l in enumerate(list_reordered)]

    #plot dN/dS
    print('plotting ds')
    OG_dict = loop_codeml(argv[1])
    dnds_list, dnds_labels = toplot_dnds(OG_dict, grouped_dict)
    list_reordered, label_list_reordered = [], []
    for k in ['busco', 'effector', 'core', 'softcore', 'accessory']:
        index = dnds_labels.index(k)
        list_reordered.append(dnds_list[index])
        label_list_reordered.append(k)
    flat_len = plot_box(list_reordered, label_list_reordered, 'dNdS_violin', True)
    for i1,i2 in itertools.combinations(range(len(list_reordered)), 2):
        print('Length mannwithney U between',
                label_list_reordered[i1], f'mean: {np.mean(list_reordered[i1])}', \
                label_list_reordered[i2], f'mean: {np.mean(list_reordered[i2])}',
                sp.stats.mannwhitneyu(list_reordered[i1], list_reordered[i2]))
    [print(f"between {label_list_reordered[i]} {np.median(l)} and all {np.median(flat_len)}: \
    {sp.stats.mannwhitneyu(l, flat_len)}") for i,l in enumerate(list_reordered)]
    #pd.DataFrame(dnds_list).to_csv('dnds_calculated.csv')
    #print(dnds_labels)

