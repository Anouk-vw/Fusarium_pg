from Bio.Phylo.PAML import codeml
from sys import argv
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

def parse_codeML(codeML_in, field):
	omega =[]
	try:
		results = codeml.read(codeML_in)
		keys =  results['pairwise'].keys()
		for kcomp in keys:
			subdict = results['pairwise'][kcomp]
			for k in keys:
				if k != kcomp:
					val = subdict[k][field]
					if val < 5:
						omega.append(val)
	except ValueError:
		#print(codeML_in)
		return False
	except KeyError:
		print(codeML_in)
		return False
	return omega#.values()

def loop(in_dir, field, end_prefix, ext):
	results_dict = {}
	for in_file in os.listdir(in_dir):
		omega_list = parse_codeML(f'{in_dir}{in_file}', field)
		if omega_list:
			OG = '_'.join(in_file.replace(ext,'').split('_')[0:end_prefix])
			results_dict[OG] = np.mean(omega_list)
	return results_dict

if __name__ == '__main__':
	#print(np.mean(parse_codeML(argv[1])))
	OG_dict_TR4 = loop(argv[1], 'omega', 2, '')
	print(np.median([float(x) for x in OG_dict_TR4.values()]))
	OG_dict_CR1TR4 = loop(argv[2], 'omega', 2, '')
	print(np.median([float(x) for x in OG_dict_CR1TR4.values()]))
	plt.boxplot([list(OG_dict_TR4.values()), list(OG_dict_CR1TR4.values())], showfliers = False)
	plt.show()
	sns.kdeplot(np.array(list(OG_dict_TR4.values())))
	sns.kdeplot(np.array(list(OG_dict_CR1TR4.values())))
	plt.show()

