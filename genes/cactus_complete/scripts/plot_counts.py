from sys import argv
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

def read_counts(counts_file):
	#read counts file from an do; echo ; wc -l; done
	count_dict = {}
	with open(counts_file) as cf:
		for line in cf:
			count_dict[line.strip()] = [int(next(cf).strip())]
	return count_dict

def add_file(count_dict, file_add):
	with open(file_add) as fa:
		for line in fa:
			count_dict[line.strip()].append(int(next(fa).strip()))
	return count_dict

def plot_dict(count_dict):
	xaxis = []
	yaxis1 = []
	yaxis2 = []
	for k,v in count_dict.items():
		xaxis.append(k)
		yaxis1.append(v[0])
		try:
			yaxis2.append(v[1])
		except:
			print(k)
			yaxis2.append(0)
	x = np.arange(len(xaxis))
	plt.bar(x-0.2, yaxis1, width=0.4, color='black')
	plt.bar(x+0.2, yaxis2, width=0.4, color='grey')
	plt.xticks(x, xaxis, rotation = 90)
	plt.ylim(14000,18700)
	plt.show()
	return

def read_gff(gff_file):
	g_lengths = []
	gene_count = 0
	with open(gff_file) as gfile:
		for line in gfile:
			if line.startswith('#'):
				continue
			if line.split('\t')[2] == 'gene':
				gene_len = int(line.split('\t')[4]) - int(line.split('\t')[3])
				g_lengths.append(gene_len)
				#if gene_len > 10000:
				#	print(line)
				gene_count += 1
			#if line.split('\t')[2] == 'gene':
			#	#g_lengths.append(gene_len)
			#	gene_count += 1
	return g_lengths, gene_count

def read_all(paths_list):
	stats = {}
	for p in paths_list:
		gname = p.split('/')[-1].split('.gff3')[0].replace('_combined', '').replace('_masked', '')
		len_list, gcount = read_gff(p)
		#print(min(len_list), len_list.index(max(len_list)))
		if np.mean(len_list) == 'NaN':
			print(len_list)
		stats[gname] = [np.mean(len_list)]
		stats[gname].append(gcount)
		#lengths.append(len_list)
	return stats

		
if __name__ == '__main__':
	#c_dict = read_counts(argv[1])
	#c_dict = add_file(c_dict, argv[2])
	#plot_dict(c_dict)
	CATgffs = pd.DataFrame.from_dict(read_all(['{}/{}'.format(argv[1], x)\
	 for x in os.listdir(argv[1]) if x.endswith('gff3') \
	 and not 'C068' in x and not 'TR4_II5' in x and not '36102' in x and not 'CR1.1' in x and not 'MINIGRAPH' in x])).T
	Fungffs = pd.DataFrame.from_dict(read_all(['{}/{}'.format(argv[2], x) \
	for x in os.listdir(argv[2]) if x.endswith('gff3_converted') \
	and not 'C068' in x and not 'TR4_II5' in x and not '36102' in x and not 'CR1.1' in x and not 'MINIGRAPH' in x])).T
	CF_df = pd.concat([CATgffs, Fungffs], axis = 1)
	CF_df.columns = ['CAT_len', 'CAT_count', 'FUN_len', 'FUN_count']
	print(CF_df.mean(axis = 0))
	CF_df[['CAT_len', 'FUN_len']].plot.bar()
	#plt.boxplot(CATgffs, positions=np.arange(len(CATgffs))-0.5, showfliers = False)
	#plt.boxplot(Fungffs, positions=np.arange(len(Fungffs)),  showfliers = False)
	#plt.show()
	CF_df[['CAT_count', 'FUN_count']].plot.bar()
	#plt.boxplot(CATgffs, positions=np.arange(len(CATgffs))-0.5, showfliers = False)
	#plt.boxplot(Fungffs, positions=np.arange(len(Fungffs)),  showfliers = False)
	plt.show()
