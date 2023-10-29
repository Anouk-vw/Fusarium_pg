
"""
gather genes per category per genome. Use bedtools to determine if within 
accessory genome.

python scripts/gene_location.py brocoli_2_noMINI/dir_step3/ gene_location/
"""

import itterate_counts as itt
from sys import argv
import numpy as np
import pandas as pd
import subprocess
from os.path import exists
import os
import re
import pybedtools as bed_tool

#def get OG classes
def add_category(counts, cutoff):
    dataframe = pd.read_csv(counts, sep = '\t', index_col=0)
    cat =[]
    for x in list((dataframe!= 0).sum(axis = 1)):
        x = int(x)
        if int(x) >= cutoff:
            cat.append('core')
        if x < cutoff and x >= cutoff-10:
            cat.append('softcore')
        if x < cutoff-10 and x>=2:
            cat.append('accessory')
        if x < 2:
            cat.append('unique')
    dataframe['category'] = cat
    return dataframe

def add_count(counts, cutoff):
    dataframe = pd.read_csv(counts, sep = '\t', index_col=0)
    count =[]
    for x in list((dataframe!= 0).sum(axis = 1)):
        if x <= 69:
            x = int(x)
            count.append(x)
        else:
            count.append(69)
    dataframe['count'] = count
    return dataframe

def OG_classes(OG_counts, num_genomes):
    """
    Read OG table and classify OGs (based on classify genes script)
    """
    tableCounts = pd.read_csv(argv[1], sep = '\t', index_col=0)
    catdf = itt.add_category(OG_counts, num_genomes)
    return catdf

#get genes per OG
def genes_OG(b_table_names, genome_name, OGs_to_list):
    """
    Read OG_names table, to return {OG: genes in genome}. Only OGs in list
    """
    #get all genes in an OG
    genes_OG = []
    with open(b_table_names) as broc_names:
        for line in broc_names:
            if line.split('\t')[0] in OGs_to_list:
                OG = line.split('\t')[0]
                genes_from_genome = [O for O in \
                    line.strip().split('\t')[1:] if genome_name in O]
                if len(genes_from_genome) > 0:
                    ##Flatten list if OG has more genes form genome
                    genes_OG.extend(list(np.concatenate([gene.split(' ') \
                    for gene in genes_from_genome])))
    return genes_OG

#group genes on class
def get_cat_genes(OG_cat, names_table, category, genome_name):
    genes = []
    OGs = OG_cat[OG_cat['count'] == category].index
    genes = genes_OG(names_table, genome_name, OGs)
    return genes

#check overlap of genes in class with adaptive genomic bedfile
def get_coordinates(genes, GFF, outfile):
    """
    return bedtools object with gff elemens occuring in adaptive regions
    """
    with open(GFF) as gff:
        generic_re = re.compile('|'.join(genes))
        matches = [line.strip().split('\t') for line in gff if bool(generic_re.search(line))]
        bed_list = [[m[i] for i in [0,3,4]] for m in matches if len(m) > 5]
        bed_object = bed_tool.BedTool.from_dataframe(pd.DataFrame(bed_list))
    return bed_object
    
def check_overlap(cat_bed, outfile):
    """
    Out: of the files in out_file, out occur in adaptive. (cat_bed = adaptive locations)
    return: [category genes in adaptive, total cat genes]
    """
    a = bed_tool.BedTool(cat_bed)
    b = bed_tool.BedTool(outfile)
    #print(len(a + b) / len(a))
    return [len(a + b) , len(a)]

def lists(info, g_core, g_acc, index, num):
    g_core.append(info[0])
    g_acc.append(info[2])
    index.append(num)
    return g_core, g_acc, index

def gene_locations(brocoli_dir, category, genome_name, outdir, gff_file\
    , acc_regions_file):
    #out_file = f'{outdir}{genome_name}_{category}_genes.bed'
    #print(out_file)
    #add_count: per count 
    #add_category: per category
    category_df = \
    add_count(f"{brocoli_dir}/table_OGs_protein_counts.txt", 69)
    #genes_per_OG = genes_OG(argv[2], 'TR4_II5')
    list_genes = get_cat_genes(category_df, \
    f"{brocoli_dir}/table_OGs_protein_names.txt", category, genome_name)
    #if not exists(out_file):
    return get_coordinates(list_genes, gff_file, out_file)

if __name__ == '__main__':
    """
    Argv[1] = brocoli dir, argv[2] = outdir, argv[3] = gff_file, argv[4] = acc_regions_file
    """
    out_d = argv[2]
    dict_locs = {}
    genomes = [x.split('.gff3')[0] for x in os.listdir(f"CAT/final_gff") if not 'Bif' in x and not 'MINI' in x]
    for cat in range(3, 69): 
        g_core = []
        g_acc = []
        index = []
        totals = 0
        adaptives = 0
        for genome in genomes:
            gff_file = f"CAT/final_gff/{genome}.gff3"
            acc_file = f"gene_location/bed_files/{genome}.bed"
            out_file = f'{out_d}{genome}_{cat}_genes.bed'
            out_bed = gene_locations(argv[1], cat, genome, out_d, gff_file, acc_file)
            adaptive, total = check_overlap(acc_file, out_bed)
            totals += total
            adaptives += adaptive
        dict_locs[cat] = adaptives/totals
    #print(dict_locs)
    np.save('fraction_adaptive.npy', dict_locs)
