from sys import argv
import os
from Bio.Phylo.PAML import codeml
import numpy as np

def parse_codeML(codeML_in):
    omega =[]
    try:
        results = codeml.read(codeML_in)
        keys =  results['pairwise'].keys()
        for kcomp in keys:
            subdict = results['pairwise'][kcomp]
            for k in keys:
                if k != kcomp:
                    #ds < 1, omega < 999
                    if float(subdict[k]['omega']) < 999:
                        omega.append(float(subdict[k]['omega']))
    except ValueError:
        #print(codeML_in)
        return False
    except KeyError:
        #print(codeML_in)
        return False
    return omega#.values()

def loop_codeml(in_dir):
    results_dict = {}
    for in_file in os.listdir(in_dir):
        omega_list = parse_codeML(f'{in_dir}{in_file}')
        if omega_list:
            OG = '_'.join(in_file.split('_')[0:2])
            results_dict[OG.replace('.txt','')] = np.mean(omega_list)
    return results_dict

def parse_deseq(deseq):
    deseq_dict = {}
    with open(deseq) as dseqf:
        for line in dseqf:
            transcript = line.split(',')[-1].replace('"', '')
            padj = line.split(',')[-2]
            l2f = line.split(',')[3]
            deseq_dict[transcript] = [l2f, padj]
    return deseq_dict

if __name__ == '__main__':
    with open('dnds_per_OG.txt', 'w+') as dnds_out:
        for k,v in loop_codeml(argv[1]).items():
            dnds_out.write(f"{k}\t{v}\n")
        
