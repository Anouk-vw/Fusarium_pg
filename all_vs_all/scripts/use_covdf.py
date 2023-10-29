"""
Read cov df and create dataframe with freq per sample.
cov_df = nucmerDF.py output
"""
from sys import argv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

##
def read_covdf(covDf):
	"""
	Get covDf
	"""
	cov_pd = pd.read_csv(covDf)
	return cov_pd

def subset(covDf, chrom, start, stop):
	"""
	From the covDF get all regions within chr,start,stop of bedfile
	"""
	#chrom_df = covDf[covDf['Chrom'] == chrom]
	cov_df_subset = covDf[(covDf['Chrom'] == chrom) &(covDf['start'] > start) & (covDf['end'] < stop)].index.tolist()
	return cov_df_subset

def subset_bed(covDf, bed):
	"""
	For bed file, get the region within chr,start,stop from the covDf
	"""
	#completeDf = pd.DataFrame(columns = covDf.columns)
	indices = []
	with open(bed) as bf:
		for line in bf:
			chrom, start, stop = line.strip().split('\t')
			temp = subset(covDf, chrom, int(start), int(stop))
			#completeDf = completeDf.append(temp)
			indices.extend(temp)
	return indices #completeDf

def sort_on(headers):
	"""
	sort meta data and obtain list of TR4/R1 IsolateCodes
	"""
	meta = pd.read_csv('/home/anouk/SPP_FOCReSequencingMetadata.csv', sep = ';')
	meta['Tr4.R1'] = meta['Tr4.R1'].fillna('R1')
	R1 = list(meta[meta['Tr4.R1'] == 'R1']['IsolateCode'])
	for s in R1:
		if s in list(headers):
			i =0
		else:
			R1.pop(R1.index(s))
	TR4 = list(meta[meta['Tr4.R1'] == 'TR4']['IsolateCode'])
	for s in TR4:
		if s in list(headers):
			i =0
		else:
			TR4.pop(TR4.index(s))
	TR4.pop(TR4.index('B2'))
	return [R1, TR4]

def sum_cols(covPd, R1_cols, TR4_cols, index):
	toplot = []
	cov_acc = covPd.iloc[index]
	cov_core = cov.drop(index)
	#sum all windows for each R1 col
	toplot.append(np.divide(list((cov_acc[R1_cols] > 0.8).sum()), len(cov_acc)))
	toplot.append(np.divide(list((cov_acc[TR4_cols] > 0.8).sum()), len(cov_acc)))
	toplot.append(np.divide(list((cov_core[R1_cols] > 0.8).sum()), len(cov_core)))
	toplot.append(np.divide(list((cov_core[TR4_cols] > 0.8).sum()), len(cov_core)))
	plt.boxplot(toplot, showfliers = False)
	plt.xticks([1,2,3,4], ['R1-acc', 'TR4-acc', 'R1-core', 'TR4-core'])
	plt.show()
	print(np.mean(toplot[0]), np.mean(toplot[1]), np.mean(toplot[2]),np.mean(toplot[3]))
	return

def plot_cov(covPd,columns, index):
	"""
	window occurs in x% of R1 isolates/x% of TR4 isolates
	get occurence for all core and accessory regions and plot
	"""
	toplot = []
	#get everything in acc index
	cov_acc = covPd.iloc[index]
	#drop everythin in acc index => core
	cov_core = cov.drop(index)
	for c in columns:
		toplot.append(cov_acc[c].tolist())
	for c in columns:
		toplot.append(cov_core[c].tolist())
	plt.boxplot(toplot, showfliers = False)
	plt.xticks([1,2,3,4], ['R1-acc', 'TR4-acc', 'R1-core', 'TR4-core'])
	plt.show()
	print(np.mean(toplot[0]), np.mean(toplot[1]), np.mean(toplot[2]),np.mean(toplot[3]))
	return

def genome_avg(covPd, R1_cols, TR4_cols, index):
	"""
	Get a 'sharedness' per genome 
	= number of windows shared per genome
	"""
	toplot = []
	cov_acc = covPd.iloc[index]
	cov_core = cov.drop(index)
	print((cov_acc[R1_cols] > 0.8).sum())
	toplot.append(np.divide(list((cov_acc[R1_cols] > 0.8).sum()), len(cov_acc)))
	toplot.append(np.divide(list((cov_acc[TR4_cols] > 0.8).sum()), len(cov_acc)))
	toplot.append(np.divide(list((cov_core[R1_cols] > 0.8).sum()), len(cov_core)))
	toplot.append(np.divide(list((cov_core[TR4_cols] > 0.8).sum()), len(cov_core)))
	plt.boxplot(toplot, showfliers = False)
	plt.xticks([1,2,3,4], ['R1-acc', 'TR4-acc', 'R1-core', 'TR4-core'])
	plt.show()
	print(np.mean(toplot[0]), np.mean(toplot[1]), np.mean(toplot[2]),np.mean(toplot[3]))
	return

if __name__ == '__main__':
	cov = read_covdf(argv[1])
	acc_i = subset_bed(cov, argv[2])
	plot_cov(cov, ['freq-R1', 'freq-TR4'], acc_i)
	R1c, TR4c = sort_on(cov.columns)
	genome_avg(cov, R1c, TR4c, acc_i)
