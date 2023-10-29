"""
usefull brocoli parsers
"""

def selective_genes_OG(b_table_names, genome_name, OGs_to_list):
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

def OGnames_dict(table_names):
    """
    return {OG: genes} from brocoli_table_names
    """
    genes_OG = {}
    with open(table_names) as broc_names:
        next(broc_names)
        for line in broc_names:
            OG = line.split('\t')[0]
            #duplicate genes sep by ' ' replace(' ', '\t') flattens list
            genes = line.replace(' ', '\t').strip().split('\t')[1:]
            #only includes genes, not empty lines (len >0)
            genes_OG[OG] = [gene for gene in genes if len(gene) > 0]
    return genes_OG
