#!/usr/bin/python3

"""
script to determine ks per duplicate type MCscanX
-> KS calculated using
1. get_sequence_pairs.py (save dup pair fasta in calc_ks/seq_pairs:
python calc_ks/get_sequence_pairs.py input_mcscan/TR4vsTR4_qcov60.blast\
    ../new_protein_files/brocoli_in/TR4_II5.new_protein.fasta \
     calc_ks/seq_pairs/
2. calc_ks.sh on fastafiles in seq_pairs. Output: codeml_out/{pair}.txt
bash calc_ks.sh seq_pairs
"""
import parse_codeml as cml
from sys import argv
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import scipy as sp
import itertools

#parse codeml
def parsecodeml(in_dir):
	codeml_dict = cml.loop(argv[1],'dS', -1,'.fasta.txt')
	return codeml_dict

#get dup type
def get_dup_type(in_file):
	type_dict = {}
	with open(in_file) as dup_file:
		for line in dup_file:
			gene, gene_type = line.strip().split('\t')
			type_dict[gene] = gene_type
	return type_dict

def link_type_ds(codeml_dict, dup_type_dict):
	kde_lists = []
	type_dict = {}
	for pair, values in codeml_dict.items():
		pair1 = pair.split('_')[0:3]
		dup_type = dup_type_dict['_'.join(pair1)]
		try:
			type_dict[dup_type].append(values)
		except KeyError:
			type_dict[dup_type] = [values]
		kde_lists.append(['_'.join(pair1), dup_type, values])
	return type_dict, kde_lists

def plot_kde(klists):
	df = pd.DataFrame(klists, columns = ['gene','duptype','dS'])
	df = df.set_index('gene')
	sns.displot(df, x ='dS', kind = 'hist', hue = 'duptype', stat='density', element = 'step',common_norm=False)
	plt.show()
	sns.violinplot(data = df, x = 'dS', y = 'duptype', cut=0)
	plt.savefig('dup_kS_violin.pdf')
	plt.show()

def plot_hist(plot_dict):
	labels = []
	boxplot_list = []
	for k,v in plot_dict.items():
		boxplot_list.append(v)
		if k == '1':
			labels.append('dispersed')
		if k == '2':
			labels.append('proximal')
		if k == '3':
			labels.append('tandem')
		if k == '4':
			labels.append('segmental')
	fig, ax = plt.subplots()
	ax.violinplot(boxplot_list, showmeans=False, showmedians=True,
        showextrema=False)

	ax.xaxis.set_tick_params(direction='out')
	ax.xaxis.set_ticks_position('bottom')
	ax.set_xticks(np.arange(1, len(labels) + 1))
	ax.set_xticklabels(labels)
	ax.set_xlim(0.25, len(labels) + 0.75)
	ax.set_xlabel('Sample name')


if __name__ == '__main__':
	ds_dict = parsecodeml(argv[1])
	dup_dict = get_dup_type(argv[2])
	type_ds_dict,kdelist = link_type_ds(ds_dict, dup_dict)
	#print(type_ds_dict.keys())
	plot_hist(type_ds_dict)
	plot_kde(kdelist)
	key_list = list(type_ds_dict.keys())
	for i1,i2 in itertools.combinations(range(len(key_list)), 2):
		k1 = key_list[i1]
		k2 = key_list[i2]
		print('dS mannwithney U between',
			k1, f'mean: {np.mean(type_ds_dict[k1])}', \
				k2, f'mean: {np.mean(type_ds_dict[k2])}',
				sp.stats.mannwhitneyu(type_ds_dict[k1], type_ds_dict[k2]))
