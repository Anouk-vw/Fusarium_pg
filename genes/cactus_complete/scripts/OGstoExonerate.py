"""
Script to determine PAV of brocoli OGs in set of isolates (exonerate)
python scripts/OGstoExonerate.py \
brocoli_2_noMINI/dir_step3/table_OGs_protein_counts.txt \
brocoli_2_noMINI/dir_step3/table_OGs_protein_names.txt \
Predicted_effectors_noMINI.txt unique_genes.txt ./
"""

from Bio import SeqIO
from sys import argv
import effector_PAV as eff_PAV
import os
import pickle

def get_OG_list(counts, names, eff, uniques):
	"""
	List all OGs with more than 50% effectors
	"""
	df = eff_PAV.fill_gene_df(counts, names, eff, uniques)
	df = df[df['category'] != 'Unique']
	percentage_effectors = df[0].astype(int) / df['gene_counts'].astype(int)
	eff_OGs = percentage_effectors[percentage_effectors >= 0.5].index
	return eff_OGs

def eff_OGs_toFasta(eff_OGs, OGs_fasta_path):
	"""
	Code to create effector fasta to run with exonerate.
	eff_OGs = list of OGs to include
	OGs_fasta_path = path where OG-fasta are stored.
					f'../OGs_fasta/prot_fasta/'
	Save longest to effectors.fasta
	
	return: dictionary to map gene_IG to OG
	"""
	#create dictionary to map gene_ID back to OG
	map_dict = {}
	#list all records to write
	records_to_save = []
	#itterate over OGs
	for eo in eff_OGs:
		lengths = []
		ids = []
		seqs = []
		records = []
		#get longest record
		for record in SeqIO.parse(f'{OGs_fasta_path}/{eo}_prot.fasta', 'fasta'):
			lengths.append(len(record.seq))
			records.append(record)
		longest = records[lengths.index(max(lengths))]
		map_dict[longest.id] = eo
		records_to_save.append(longest)
	#with open("effectors.fasta", "w") as output_handle:
	#	SeqIO.write(records_to_save, output_handle, "fasta")
	with open('map_dict.pkl', 'wb') as f:
		pickle.dump(map_dict, f)
	return map_dict

if __name__ == '__main__':
	#step 1
	broc_counts = argv[1]
	broc_names = argv[2]
	eff_txt = argv[3]
	uniques_txt = argv[4]
	fasta_path = argv[5]
	OGs_list = get_OG_list(broc_counts, broc_names, eff_txt, uniques_txt)
	if os.path.exists('effectors.fasta'):
		eff_OGs_toFasta(OGs_list, fasta_path)
	#step 2 (run exonerate inbetween)
