"""
example:
python nucmerDF.py nuc_synt/sliding_windows/ CR1.1 synt_sliding.bed windows CR1.1_covdf.csv

Read nucmer alignments in bed format (windows) from Snakefile 
create dataframe with coverage per genome per window

argv[1] = directory to read
argv[2] = genome to use as reference
argv[3] = extension of files to read
argv[4] = sperator
argv[4] = location to save dataframe

Dataframe can be used to plot with circos (genome cirocs, acc circos) 
ToDo: Create these as scripts as well
"""

import pandas as pd
import os
from sys import argv

#functions
def init_df(indir, genome_id,covcutoff, ext, sep):
    """
    Read all windows files as output by nucmer snakemake
    input: 
    genome = name of 'reference' genome
    ext = ending of the files to read ('sliding.bed')
    
    output:
    dataframe with windows of reference, coverage per genome per window
    """
    header_inc = []
    cov_dfs = []
    headers = []
    inputFiles = [x for x in os.listdir(indir) if x.startswith('{}{}'.format(genome_id, sep)) and x.endswith(ext)]
    print(len(inputFiles))
    for f in inputFiles:
        #print(f)
        df = pd.read_csv('{}{}'.format(indir, f), sep = '\t', header = None)
        if len(cov_dfs) > 0:
            cov_dfs.append(df[6])
            headers.append(f.split(sep)[1].split(ext)[0])
            header_inc.append(f)
        else:
            header_inc.append(f)
            cov_dfs.append(df)
            headers.extend(['Chrom', 'start', 'end', 'counts', 'bp', 'total', f.split(sep)[1].split(ext)[0]])
    cov_df = pd.concat(cov_dfs, axis = 1)
    cov_df = pd.concat((cov_df, (cov_df.iloc[:,6:,] >= covcutoff).sum(axis = 1)), axis = 1) #(cov_df.iloc[:,6:-1].mean(axis = 1))), axis = 1)
    headers.append('count')
    cov_df.columns = headers
    print(len(header_inc))
    return cov_df

def read_meta(meta):
    """
    read meta and extract R1, TR4, R2 isolate names.
    '/home/anouk/SPP_FOCReSequencingMetadata.csv'
    """
    m = pd.read_csv(meta, sep = ',')
    m['Tr4.R1'] = m['Tr4.R1'].fillna('R1')
    R1_names = list(m[m['Tr4.R1'] == 'R1']['IsolateCode'])
    R2_names = list(m[m['Tr4.R1'] == 'R2']['IsolateCode'])
    TR4_names = list(m[m['Tr4.R1'] == 'TR4']['IsolateCode'])
    return [R1_names, R2_names, TR4_names]

def remove_missing(race_names, cov_df):
    """
    For each name from meta data, check if it is present in df cols
    """
    remove = []
    for x in race_names:
        #print(x)
        if x not in cov_df.columns:
            remove.append(race_names.index(x))
    for index in sorted(remove, reverse=True):
        del race_names[index]
    return race_names

def add_counts(cov_df, covcutoff, names):
    """
    add column with counts and synteny. Add frequency column for races.
    cov_df = dataframe with coverage from init_df
    covcutoff = cutoff of nucleotides covered to call present
    names = list of lists with isolate names. [[R1], [R2], [TR4]]
    """
    headers = list(cov_df.columns)
    cov_df = pd.concat((cov_df, (cov_df.iloc[:,6:,] >= covcutoff).sum(axis = 1)), axis = 1) #(cov_df.iloc[:,6:-1].mean(axis = 1))), axis = 1)
    headers.append('count')
    cov_df = pd.concat((cov_df, (cov_df.iloc[:,6:-1]).mean(axis = 1)), axis = 1) #(cov_df.iloc[:,6:-1].mean(axis = 1))), axis = 1)
    headers.append('synt ratio')
    for race, n in zip(['R1', 'R2', 'TR4'], names):
        h = 'freq-{}'.format(race)
        n = remove_missing(n, cov_df)
        R_df = cov_df[n]
        count = (R_df >= covcutoff).sum(axis = 1)
        cov_df = pd.concat((cov_df, count/len(n)), axis = 1) #(cov_df.iloc[:,6:-1].mean(axis = 1))), axis = 1)
        headers.append(h)
    cov_df.columns = headers
    return cov_df

def create_df(indir, genome, ext, sep, cutoff, meta):
    """
    Create final df with coverage/sample/window and counts
    genome = genome name as reference
    cutoff = coverage to call window present
    meta = location of meta data
    """
    cov_df = init_df(indir, genome, cutoff, ext, sep)
    #iso_names = read_meta(meta)
    #cov_df = add_counts(cov_df, cutoff, iso_names)
    return cov_df

def save_df(indir, genome, ext, sep, cutoff, output):
    print(indir, genome , cutoff, ext, sep, output)
    df = init_df(indir, genome , cutoff, ext, sep)
    df.to_csv(output)

if __name__ == '__main__':
    df = create_df(argv[1], argv[2], argv[3], argv[4], 0.8, '/home/anouk/SPP_FOCReSequencingMetadata.csv')
    df.to_csv(argv[5])
