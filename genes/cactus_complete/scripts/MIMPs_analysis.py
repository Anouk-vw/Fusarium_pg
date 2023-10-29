"""
Get effectors close to MIMP elements
"""
import sys
from sys import argv
sys.path.append('/home/anouk/anouk2/pangenome/cactus_complete/genes/scripts')
import effector_PAV as eff_PAV
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def get_OG_list(counts, names, eff, uniques,buscos):
	"""
	List all OGs with more than 50% effectors
	"""
	df = eff_PAV.fill_gene_df(counts, names, eff, uniques, buscos)
	df = df[df['category'] != 'Unique']
	percentage_effectors = df['effectors'].astype(int) / df['gene_counts'].astype(int)
	eff_OGs = percentage_effectors[percentage_effectors >= 0.5].index
	return eff_OGs

def create_OG_dict(b_table_names, OG_list):
    """
    Read OG_names table, to return {OG: genes in genome}. Only OGs in list
    """
    #get all genes in an OG
    genes_OG = {}
    with open(b_table_names) as broc_names:
        for line in broc_names:
            OG = line.split('\t')[0]
            if OG in OG_list:
                genes_OG[OG] = line.split('\t')[1:]
            """
            genes_from_genome = [O for O in \
                line.strip().split('\t')[1:] if genome_name in O]
            if len(genes_from_genome) > 0:
                ##Flatten list if OG has more genes form genome
                genes_OG.extend(list(np.concatenate([gene.split(' ') \
                for gene in genes_from_genome])))
            """
    return genes_OG

def parse_list_file(input_file):
    item_list = []
    with open(input_file) as list_file:
        for line in list_file:
            item_list.append(line.strip())
    return item_list

def gene_ID_vs_OG(OG_genes_list, gene_IDs):
    """
    check if gene_IDs occurs in OG
    """
    return bool(set(gene_IDs) & set(OG_genes_list))

def number_gene_ID_vs_OG(OG_genes_list, gene_IDs):
    """
    check if gene_IDs occurs in OG
    """
    return len(set(gene_IDs) & set(OG_genes_list))

def genes_in_OGs(counts_file, names_file, eff_file, unique_file, gene_ID_list, buscos):
    """
    Check per OG is gene from gene list occurs
    Broccoli output, effectorp file and unique genes
    """
    OGs_with_gene = []
    eff_list = get_OG_list(counts_file, names_file, eff_file, unique_file, buscos)
    OG_dict = create_OG_dict(names_file, eff_list)
    for OG_id, OG_genes in OG_dict.items():
        if gene_ID_vs_OG(OG_genes, gene_ID_list):
            OGs_with_gene.append(OG_id)
        else:
            continue
    return OGs_with_gene

def make_clustermap(counts, phyl_order, OG_index):
    table = eff_PAV.get_OG_cat(pd.read_csv(counts, sep = '\t', index_col=0), 69)
    table.columns = [x.replace('.new_protein.fasta', '') for x in table.columns]
    table_cluster = (table[table['category'] != 'core'] !=0).iloc[:,:69]
    ordered = table_cluster.reindex(phyl_order, axis = 1)
    subset_table = ordered.reindex(OG_index).fillna(False)
    print(subset_table)
    sns.clustermap(subset_table, figsize =(25,10))
    plt.show()
    return
    
if __name__=='__main__':
    broc_dir = argv[1]
    effector_file = argv[2]
    unique_file = argv[3]
    gene_list = parse_list_file(argv[4])
    busco = '/home/anouk/anouk2/pangenome/cactus_complete/genes/new_protein_files/busco_annotation/all_complete_busco_genes.txt'
    Ogs = genes_in_OGs(f'{broc_dir}/table_OGs_protein_counts.txt',\
                f'{broc_dir}/table_OGs_protein_names.txt',\
                effector_file,
                unique_file,
                gene_list)
    phylogeny = ['36102','BRIP62280','BRIP63632','Leb1-2C','BRIP65062',
    'M3','M2','M4','M5','S1B8','Race4','Pak1.1A','TR4_II5','II_5','M1',
    'Eden','ISR1','ISR5','JV11','JV14','Laos','Vietnam','Myanmar',
    'Col17','Col2','Col4','Indo8','Phi2-6C', '9_KT06-A1', 'NRRL_36107', 
    'NRRL_36103', 'Bif_04', 'CR1.1', 'FocST4-98','FocP1', 'NRRL_36101', 
    'NRRL_36112', 'NRRL_36110', 'C058', 'C187','C082', 'C176', 'C177', 
    'C081', 'C192', 'C135', 'C090', 'NRRL_36117','NRRL_36108', 'P41b', 
    'NRRL_36113', 'NRRL_36118', '19_KB07-B','56_JB09-A2', 'NRRL_36120', 
    'NRRL_36115', 'Mal43', 'NRRL_36116','15_KTG06-B2', 'Race1', 'P26c', 
    'P20a', 'Phi6.6a', 'Indo110']
    make_clustermap(f'{broc_dir}/table_OGs_protein_counts.txt',\
    phylogeny, Ogs)
