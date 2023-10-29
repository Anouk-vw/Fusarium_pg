"""
Script to get single genes (gene not in any cluster) from brocoli output
Input: 
argv[1] = output table. brocoli/dir_step3/table_names
argv[2] = directory. gene_predictions/gff

example usage:
python singletons_brocCol.py dir_step3/table_OGs_protein_names.txt \
./CAT/consensus_gene_set/ > singletons_counts.txt 
"""
from sys import argv
import re
from Bio import SeqIO
import os

def genes_OG(inpath, genome_name, ext):
	""". 
	Return: gene names in OG (for genome)
	"""
	names_in_OG = []
	#read table_names form brocoli
	with open(argv[1]) as brocfile:
		for line in brocfile:
			#get header
			if line.startswith('#'):
				header_list = line.strip().split('\t')
				#get index of gene_name
				hi = header_list.index(f'{genome_name}{ext}')
			#get all gene_names at index 'hi'
			gene_names = line.split('\t')[hi]
			#if gene names has string, append. Else, col[hi] is empty.
			if len(gene_names) > 1:
				names_in_OG.extend(gene_names.strip().split(' '))
	return names_in_OG

def transcript_total(inpath, genome_name):
	"""
	Return all transcripts in gff3
	"""
	all_transcripts = []
	#read gff
	with open(f"{inpath}/{genome_name}.gff3") as genes:
		for line in genes:
			#if not header and line contains transcript info
			if not line.startswith('#') and line.split('\t')[2] == 'transcript':
				#find transcript id and return id
				transcript_id = re.findall('transcript_id=(.*?);', line)
				#extend all ids
				all_transcripts.extend(transcript_id)
	return all_transcripts

def genes_total(inpath, genome_name, ext):
	"""
	Return all prot coding genes in fasta (1 transcript per gene)
	fasta created by 'scripts/gff_to_protein.py'
	"""
	all_genes = []
	genome_fasta = f'{inpath}{genome_name}{ext}'
	#read gff
	for record in SeqIO.parse(genome_fasta, 'fasta'):
		all_genes.append(record.id)
	return all_genes

def save_single_names(all_names, OGs, out):
	"""
	return: file with all transcript names not in OG.{genome}_singletons
	all_names = list of all  protein coding transcripts (fed to brocoli)
	"""
	#gene not in OG = singleton. in all names not in OG names
	singleton_names = set(all_names) - set(OGs)
	#write to file
	with open(f"{out}_singletons", 'w+') as outfile:
		for item in singleton_names:
			outfile.write(f"{item}\n")
	return

if __name__ == '__main__':
	inpath_OG = argv[1]
	inpath_total = argv[2]
	#extension is ending of protein fastas
	extension = argv[3]
	genome_name = [x for x in os.listdir(argv[2]) if \
			x.endswith(extension) and not x.startswith('_')]
	for g in genome_name:
		genome_name = g.replace(extension, '')
		OG_genes = genes_OG(inpath_OG, genome_name, extension)
		total_genes = genes_total(inpath_total, genome_name, extension)
		save_single_names(total_genes, OG_genes, genome_name)
		print(f"{genome_name}\t{len(set(total_genes)) - len(set(OG_genes))}")
