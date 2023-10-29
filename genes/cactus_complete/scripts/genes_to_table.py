"""
python genes_to_table.py ../../../TR4_loc_categories.bed ../kallisto/quantification/DEseq_results8dpi_TR4II5.csv gene_age.txt 
python genes_to_table.py ../../../TR4_loc_categories.bed ../kallisto/quantification/DEseq_results8dpi_TR4II5.csv gene_age_atleast1.txt dnds_per_OG.txt ../MCscanX_dups/input_mcscan/TR4vsTR4_qcov60.gene_type 
"""

import sys
sys.path.append('/home/anouk/anouk2/pangenome/cactus_complete/genes/')
import scripts.MIMPs_analysis as MIMP
import scripts.effector_PAV as eff_PAV
import pandas as pd
from sys import argv

#def read_input(counts_table, names_table, cutoff, ext): 
#    """
##    cutoff = Int. Occurences to consider core 
##    ext = string to replace. Can be ''
#    """
##    table = eff_PAV.get_OG_cat(pd.read_csv(counts_table, sep = '\t', index_col=0), cutoff)
#    table.columns = [x.replace(f'{ext}', '') for x in table.columns]
#    return table

def link_OG_genes_OrthoFinder():
    OG_gene_names = {}
    with open('Results_Feb23_1/Orthogroups/Orthogroups.tsv') as OG_counts_file:
        for line_num, line in enumerate(OG_counts_file):
            if line_num == 0:
                header = line.split('\t')
            #if line_num > 0 and line_num < 10:
            else:
                OG_gene_names[line.split('\t')[0]] =  \
                {header[i+1]:x for i,x in enumerate(line.split('\t')[1:])}
    return OG_gene_names

def map_OG_TR4(names_file, genome):
    OG2TR4 = {}
    with open(names_file) as names:
        for line in names:
            if not line.startswith('#'):
                for x in line.split('\t'):
                    if 'OG' in x:
                        OG = x
                    if genome in x:
                        TR4_gene = x
                OG2TR4[OG] = TR4_gene
    return OG2TR4

def get_geneloc_dict(gp):
    """
    dict {gene:loc} loc as chr-start-stop
    gp_file = path were .gp files are stores
    """
    gp_dict = {}
    with open(gp) as gp_file:
        for line in gp_file:
            line = line.strip()
            gene = line.split('\t')[0]
            chrom = line.split('\t')[1]
            start = line.split('\t')[3]
            stop = line.split('\t')[4]
            gp_dict[gene] = [chrom, start, stop] #f'{chrom}-{start}-{stop}'
    return gp_dict

def parse_info_file(info_brocoli):
    info_dict = {}
    with open(info_brocoli) as info_file:
        for line in info_file:
            if not line.startswith('#'):
                OG = line.strip().split(',')[0]
                info = line.strip().split(',')[1:]
                info_dict[OG] = info
    return info_dict

def get_loc_cat(bed_file):
    cat_loc_dict = {}
    with open(bed_file) as bfile:
        for line in bfile:
            chrom, start, stop, cat = line.strip().split('\t')
            try:
                cat_loc_dict[chrom][f'{start}-{stop}'] = cat
            except KeyError:
                cat_loc_dict[chrom] = {}
                cat_loc_dict[chrom][f'{start}-{stop}'] = cat
    return cat_loc_dict

def check_loc(location, location_dict):
    #if loc in location file and loc < stop
    #else: cat is core
    chrom, start, stop = location
    for k,v in location_dict[chrom].items():
        coord_cat_start, coord_cat_stop = k.split('-')
        if int(stop) > int(coord_cat_start) and int(start) < int(coord_cat_stop):
            #print(chrom, start, stop, coord_cat_start, coord_cat_stop)
            return v
        else:
            return 'core'

def parse_deseq(deseq):
    deseq_dict = {}
    with open(deseq) as dseqf:
        for line in dseqf:
            line = line.strip()
            transcript = line.split(',')[-1].replace('"', '')
            padj = line.split(',')[-2]
            l2f = line.split(',')[3]
            deseq_dict[transcript] = [l2f, padj]
    return deseq_dict

def parse_dnds(dnds_file):
    dnds_dict = {}
    with open(dnds_file) as dnds_f:
        for line in dnds_f:
            OG, dnds_mean, dnds_median = line.strip().split('\t')
            dnds_dict[OG] = dnds_median
    return dnds_dict

def parse_gene_age(age_file):
    age_dict = {}
    with open(age_file) as age_f:
        for line in age_f:
            transcript, phylum, rank = line.strip().split('\t')
            age_dict[transcript] = [phylum, rank]
    return age_dict

def parse_dup_type(dup_file):
    #../MCscanX_dups/input_mcscan/TR4_II5vsTR4_II5_qcov50.gene_type 
    dup_dict = {}
    with open(dup_file) as dups:
        for line in dups:
            gene, dup_type = line.strip().split('\t')
            dup_dict[gene] = dup_type
    return dup_dict

def fill_gene_table(gp_path, genome, broc_names, info_brocoli, loc_dict, deseq_d, age_dict, dnds_dict, dup_dict):
    init_dict = {}
    #get {OG:[genes]}
    OG_dict = map_OG_TR4(broc_names, genome)
    #get {gene:loc}
    gene_location = get_geneloc_dict(gp_path)
    #init_dict = {k:[v] for k,v in gene_location.items()}
    #for all info, append info to appropriate gene
    info_dict = parse_info_file(info_brocoli)
    for OG, info_OG  in info_dict.items():
        try:
            dnds_val = dnds_dict[OG]
        except KeyError:
            dnds_val = 'NA'
        if OG in OG_dict:
            genes = OG_dict[OG].split(' ')
            for gene in genes:
                info = []
                info.extend(info_OG)
                #get location
                info.extend(gene_location[gene])
                #categorize location
                info.append(check_loc(gene_location[gene], loc_dict))
                #get deseq
                try:
                    info.extend(deseq_d[gene])
                except KeyError:
                    info.extend(['NA', 'NA'])
                    #print(gene)
                #add age
                try:
                    info.extend(age_dict[gene])
                except KeyError:
                    info.extend(['NA', 'NA'])
                #add dnds
                info.append(dnds_val)
                #add dup type
                info.append(dup_dict[gene])
                
                init_dict[gene] = info
        else:
            #OG is unique
            if OG.startswith(genome):
                gene = OG
                info = []
                info.extend(info_OG)
                #get location
                info.extend(gene_location[gene])
                #categorize location
                info.append(check_loc(gene_location[gene], loc_dict))
                #get deseq
                try:
                    info.extend(deseq_d[gene])
                except KeyError:
                    info.extend(['NA', 'NA'])
                #add age
                try:
                    info.extend(age_dict[gene])
                except KeyError:
                    info.extend(['NA', 'NA'])
                #add dnds
                info.append(dnds_val)
                #add dup type
                info.append(dup_dict[gene])
                
                init_dict[gene] = info
            #print('problem with OG:', OG)
    return init_dict
    

if __name__ == '__main__':
    brocoli_path  = '/home/anouk/anouk2/pangenome/cactus_complete/genes/brocoli_results'
    counts = f'{brocoli_path}/dir_step3/table_OGs_protein_counts.txt'
    names = f'{brocoli_path}/dir_step3/table_OGs_protein_names.txt'
    eff = '/home/anouk/anouk2/pangenome/cactus_complete/genes/Predicted_effectors_noMINI.txt'
    uniques = '/home/anouk/anouk2/pangenome/cactus_complete/genes/brocoli_results/unique_genes.txt'
    busco_genes = '/home/anouk/anouk2/pangenome/cactus_complete/genes/new_protein_files/busco_annotation/all_complete_busco_genes.txt'
    gp_path = "/home/anouk/anouk2/pangenome/cactus_complete/CAT/consensus_gene_set/TR4_II5.gp"
    #read orthogroups brocoli
    #broc_table = read_input(counts, names, 69, '.new_protein.fasta')
    #read categories brocoli
    #For now ignore Unique
    #info_df = pd.read_csv('info_broccoli.csv', index_col = 0)
    #info_df = info_df[info_df['category'] != 'Unique']
    
    #init_df
    location_cat_dict = get_loc_cat(argv[1])
    deseq = parse_deseq(argv[2])
    age = parse_gene_age(argv[3])
    dnds = parse_dnds(argv[4])
    dups = parse_dup_type(argv[5])
    
    init_gene_dict = fill_gene_table(gp_path , 'TR4_II5', names, '../info_broccoli.csv', location_cat_dict, deseq, age, dnds, dups)
    
    #print({k:v for k,v in init_gene_dict.items() if len(v) != 16})
    
    pd.DataFrame.from_dict(init_gene_dict).T.to_csv('gene_info_table.csv')
    
    #print(check_loc(['TR4_II5_Chr1', 0, 100], get_loc_dict(argv[1])))
