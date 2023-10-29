from sys import argv
from Bio import SeqIO

def get_sequence(OGs, prot_dir):
    longest = []
    for OG in OGs:
        lengths = []
        records = []
        for record in SeqIO.parse(f'{prot_dir}/{OG}_prot.fasta', 'fasta'):
            lengths.append(len(record.seq))
            records.append(record)
        longest.append(records[lengths.index(max(lengths))])
    [print(f'>{l.id}\n{l.seq}') for l in longest]
    return 
    
def genes_from_OG(genome, names_df, gene_list):
    TR4_genes_OG = names_df[names_df[f"{genome}.protein.consensus.fasta"].notna()][f"{genome}.protein.consensus.fasta"]
    print(TR4_genes_OG.reindex(gene_list))
    acc_genes_TR4 = TR4_genes_OG.reindex(gene_list)[TR4_genes_OG.reindex(gene_list).notna()]
    return acc_genes_TR4
    
if __name__ == '__main__':
    OGs = argv[3:]
    prot_loc = argv[2]
    if argv[1] == 'get_sequence':
        print(get_sequence(OGs, prot_loc))
    if argv[1] == 'gene_ids':
        names_table = argv[3]
        print(genes_OG(names_table, OGs))
    
