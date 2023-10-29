from sys import argv
import pandas as pd
import random
import numpy as np
import itterate_counts as itt
import seaborn as sns
import matplotlib.pyplot as plt

"""
argv[1] = table_counts
argv[2] = table_names
argv[3] = effector_gene_names
(effector gene_names:
grep -v 'Non' effectorp/* | grep -P '_T[0-9]*' | cut -f 1 | cut -d':' -f 2 >> Predicted_effectors.txt)
"""

def get_OG_cat(broc_df, cutoff):
    """
    Return dataframe with OG-genomes PAV-cate
    """
    dataframe = itt.add_category(broc_df, cutoff)
    return dataframe

def genes_OG(b_table_names):
    #get all genes in an OG
    genes_OG = {}
    with open(b_table_names) as broc_names:
        for line in broc_names:
            if line.startswith('#'):
                continue
            else:
                OG = line.split('\t')[0]
                ##This returns a list if one genome has more than one gene in teh cluster. Flatten list
                genes_OG[OG] = list(np.concatenate([gene.split(' ') for gene in \
                        line.strip().split('\t')[1:] if len(gene) > 1]))
    return genes_OG

def effector_list(effector_names):
    eff_list = []
    with open(effector_names) as eff_file:
        for line in eff_file:
            eff_list.append(line.strip())
    return set(eff_list)

def is_gene_effector(gene_name, effector_list):
    return gene_name in effector_list

def OG_effectors(OG_genes_dict, effector_list):
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

def core_acc_dict(OGs, dataframe):
    ca_dict = {}
    core = dataframe.index[dataframe['category'] == 'core']
    acc = dataframe.index[dataframe['category'] == 'accessory']
    sc = dataframe.index[dataframe['category'] == 'softcore']
    ca_dict['core'] = list(np.concatenate([OGs[OG] for OG in core]))
    ca_dict['acc'] = list(np.concatenate([OGs[OG] for OG in acc]))
    ca_dict['softcore'] = list(np.concatenate([OGs[OG] for OG in sc]))
    return ca_dict

def eff_count_from_dict(CA_dict,  effector_list):
    core_eff = []
    acc_eff = []
    sc_eff = []
    for k,v in CA_dict.items():
        if k == 'core':
            for gene in v:
                core_eff.append(is_gene_effector(gene, effector_list))
        if k == 'softcore':
            for gene in v:
                sc_eff.append(is_gene_effector(gene, effector_list))
        if k == 'acc':
            for gene in v:
                acc_eff.append(is_gene_effector(gene, effector_list))
    return core_eff, sc_eff, acc_eff

def eff_occurence(CA_dict,  effector_list):
    core_eff = 0
    acc_eff = 0
    for e in effector_list:
        if e in set(CA_dict['core']):
            core_eff += 1
        if e in set(CA_dict['acc']):
            acc_eff +=1
    return core_eff, acc_eff

def get_unique(unique, effector_list):
    eff_in_u = []
    u_eff = {}
    with open(unique) as u:
        for line in u:
            u_eff[line.strip()] = is_gene_effector(line.strip(), effector_list)
            #eff_in_u.append(is_gene_effector(line.strip(), effector_list))
    return u_eff

if __name__ == '__main__':
    tableCounts = pd.read_csv(argv[1], sep = '\t', index_col=0)
    catdf = get_OG_cat(tableCounts, 69)
    OG_dict = genes_OG(argv[2])
    effectors = effector_list(argv[3])
    OG_eff_dict = OG_effectors(OG_dict, effectors)
    OG_eff_dict = {k:1 if any(v) else 0 for k,v in OG_eff_dict.items()}
    ca = core_acc_dict(OG_dict, tableCounts)
    #eff_OGs = len([k for k,v in OG_eff_dict.items() if any(v)])
    u_dict = get_unique(argv[4], effectors)
    udf = pd.DataFrame.from_dict(u_dict, orient = 'index')
    odf = pd.DataFrame.from_dict(OG_eff_dict, orient = 'index')
    concat = pd.concat([udf, odf])
    comb_df = pd.concat([tableCounts, concat], axis = 1)
    comb_df['category'] = comb_df['category'].replace(np.nan, 'Unique')
    print(comb_df)
    eff_core = comb_df[['category', 0]].value_counts(normalize = True)
    print(comb_df)
    eff_core.plot.bar()
    plt.show()
    """
    c_eff, soft_eff, a_eff = eff_count_from_dict(ca, effectors)
    effu = get_unique(argv[4], effectors)
    eff_c_u = (np.sum(effu) / len(effectors)) * 100
    eff_c_c = (np.sum(c_eff)/len(effectors))*100
    eff_c_a = (np.sum(a_eff)/len(effectors))*100
    eff_c_sc = (np.sum(soft_eff)/len(effectors))*100
    print(len(effectors), np.sum(effu) + np.sum(c_eff) + np.sum(a_eff) + np.sum(soft_eff))
    labels = ['core', 'softcore','acc','unique']
    values  = [eff_c_c, eff_c_sc, eff_c_a, eff_c_u]
    fig1, ax1 = plt.subplots()
    ax1.pie(values, labels=labels, autopct='%1.1f%%',
        shadow=True, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.show()
    #c_eff, a_eff = eff_occurence(ca, effectors)
    #print(c_eff)
    """
