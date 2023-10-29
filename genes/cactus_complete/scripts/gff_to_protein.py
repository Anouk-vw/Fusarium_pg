"""
script to create a protein file with one representative transcript per gene

"""
import pprint
from BCBio import GFF
from Bio import SeqIO
from sys import argv
import re
import os

def gene_id_transcript_dict(gff_file):
    """
    Link gene_id to transcript_ids. Collect isoforms/gene
    input: gff_file
    output: {gene_id:[transcripts]}
    """
    gene_transcript_dict = {}
    in_file = gff_file
    with open(in_file) as gff:
        for line in gff:
            """
            if not '#' in line and line.split('\t')[2] == 'gene':
                #search for 'gene_id=.*;' and get gene_id
                result = re.findall('gene_id=(\w*);', line)
                gene_id = result[0]
            """
            if not '#' in line and line.split('\t')[2] == 'transcript':
                result = re.findall('gene_id=([\w\.\-]*);', line)
                try:
                    gene_id = result[0]
                except:
                    print(result, gff_file, line)
                result = re.findall('transcript_id=([\w\.\-]*);', line)
                transcript_id = result[0]
                try:
                    gene_transcript_dict[gene_id].append(transcript_id)
                except:
                    gene_transcript_dict[gene_id] = [transcript_id]
    return gene_transcript_dict

def get_lengths(fasta_file):
    """
    fasta_file = cat out
    output: {transcript_id:length}
    """
    id_len_dict = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        id_len_dict[record.id] = len(record.seq)
    return id_len_dict

def get_longest(trans_len_dict, transcripts_list):
    lengths = []
    for T in transcripts_list:
        try:
            lengths.append(int(trans_len_dict[T]))
        except KeyError:
            #some T's are already ommited from the consensus apperently
            print(T)
    #lengths = [int(trans_len_dict[T]) for T in transcripts_list]
    longest_id = lengths.index(max(lengths))
    longest_T = transcripts_list[longest_id]
    return longest_T

def get_transcripts_to_include(gene_transcript_dict, transcript_len):
    """
    Select longest transcript
    return: {gene_id:longest_transcript}
    """
    final_transcript_dict = gene_transcript_dict
    for gene, transcripts in gene_transcript_dict.items():
        #more than 1 transcript per gene
        if len(transcripts) > 1:
            #save longest transcript
            final_transcript_dict[gene] = get_longest(transcript_len, transcripts)
        else:
            final_transcript_dict[gene] = transcripts[0]
    return final_transcript_dict

def create_fasta(trans_toinclude_dict, original_fasta, out_fasta):
    to_save = []
    gene_ids, trans_ids = zip(*trans_toinclude_dict.items())
    for record in SeqIO.parse(original_fasta, 'fasta'):
        if record.id in trans_ids:
            to_save.append(record)
    """print(len(to_save))
    print(len(gene_ids))
    print(set(trans_ids) - set([x.id for x in to_save]))"""
    with open(out_fasta, "w") as output_handle:
        SeqIO.write(to_save, output_handle, "fasta")
    return

def gff_to_fasta(gff, fasta, out):
    gene_dict = gene_id_transcript_dict(gff)
    lengths = get_lengths(fasta)
    final_dict = get_transcripts_to_include(gene_dict, lengths)
    create_fasta(final_dict, fasta, out)
    return

if __name__ == '__main__':
    indir = argv[1] #CAT/consensus_gene_set
    outdir = argv[2]
    #get all items in CAT_directory
    genomes = [x.replace('.gff3', '') for x in os.listdir(indir) if '.gff3' in x]
    #itterate over items
    for genome in genomes:
        print(f'working on genome {genome}')
        prefix = f'{indir}{genome}'
        #don't create new file if genome.new_protein.fasta, exists
        if os.path.exists(f'{outdir}/{genome}.new.fasta'):
            print(f'{genome} exists')
            continue
        else:
        #create fasta from gff and write to x/x.new_protein.fasta
            gff_to_fasta(f'{prefix}.gff3', \
                        f'{prefix}.consensus.fasta', \
                        f'{outdir}/{genome}.new.fasta')
