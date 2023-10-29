"""
Analyse nucmer output from snakefile
input: bed files with first name as reference
output: 
1. Heatmap of co-occuring adaptive regions between genomes
2. Size of acc/core content per genome
...

python analyse_nucmer.py cov_df_out cov_df_sliding Plots

argv[1] = path to directory with cov_dfs
argv[2] = path to directory with sliding windows cov_df
argv[3] = prefix for output figures
"""

from sys import argv
import pandas as pd 
import os
import matplotlib.pyplot as plt
import seaborn as sns

#Read dataframes
def read_df(dataframe_path):
    """
    read dataframe in path
    """
    df = pd.read_csv(dataframe_path)
    return df

def add_core_dict(df, cutoffs):
    """
    Obtain genome size for provided dataframe
    input: df = dataframe (row: ref-windows, col: genomes)
           cutoffs = cutoff to consider core
    output: list of sizes. core, softcore, accessory, unique
    """
    core = df[df['count'] >= cutoffs]['total'].sum()
    softcore = df[df['count'].between(cutoffs*0.8, cutoffs)]['total'].sum()
    accessory = df[(df['count'].between(1, cutoffs*0.8))]['total'].sum()
    unique = df[(df['count'] == 0)]['total'].sum()
    #current_dict[genome_name] = [core, softcore, accessory, unique]
    return [core, softcore, accessory, unique]

def fill_sizes_dict(genome_list, path_to_covdf, cutoff):
    """
    Add genome size to dict for all genomes in genome_list
    output: dict{genome:[sizes]}
    """
    #genome_list = [x.split('_covdf')[0] for x in os.listdir(path_to_covdf)]
    #cutoff= len(genome_list) - 2
    #print(cutoff)
    genomes_dict = {}
    for genome_name in genome_list:
        genomevs_df = read_df(f'{path_to_covdf}/{genome_name}_covdf.csv')
        #genomes_dict = fill_core_dict(genomevs_df, genomes_dict, genome_name, cutoff)
        genomes_dict[genome_name] = add_core_dict(genomevs_df, cutoff)
    return genomes_dict
    
def order_meta(meta, sizes_dict, order_on, return_col):
    """
    use metadata to order the species. Only keep meta_labels in dict
    input:
    order_on = columns in meta dat to order on
    return_col = column in meta data to return
    """
    meta['Tr4.R1'] = meta['Tr4.R1'].fillna('R1')
    ordered_species = list(meta.sort_values([order_on], axis =0)[return_col])
    for s in ordered_species:
        try: 
            sizes_dict[s]
        except:
            ordered_species.pop(ordered_species.index(s))
    return ordered_species

def plot_sizes(genomes_dict, out, ordered_species=None):
    """
    return plot of core/acc sizes per genome
    input: genomes_dict = dict{genome_name:[sizes]}
    ordered_species: optional. Provide list to order the x axis
    """
    if ordered_species == None:
        ordered_species = genomes_dict.keys()
        print(ordered_species)
    plt.subplots(figsize=(20,5))
    plt.bar(ordered_species, [v[0] for v in list(genomes_dict.get(i) for i in ordered_species)], color = '#07738F')
    plt.bar(ordered_species, [v[1] for v in list(genomes_dict.get(i) for i in ordered_species)], bottom = [v[0] for v in list(genomes_dict.get(i) for i in ordered_species)], color = '#16B3DB')
    plt.bar(ordered_species, [v[2] for v in list(genomes_dict.get(i) for i in ordered_species)], bottom = [v[1] + v[0] for v in list(genomes_dict.get(i) for i in ordered_species)], color = '#DB004C')
    plt.bar(ordered_species, [v[3] for v in list(genomes_dict.get(i) for i in ordered_species)], bottom = [v[1] + v[0] +v[2] for v in list(genomes_dict.get(i) for i in ordered_species)], color = '#DBCB16')
    plt.legend(['core', 'softcore', 'accessory', 'unique'])
    plt.xticks(rotation = 90)
    #plt.show()
    plt.savefig(f'{out}_sizes.pdf')

def get_shared_per_genome(genome, above, below, path_to_cov):
    """
    Get mean shared acc_regions per genome for a single ref 
    input: 
    genome = str(reference name)
    above, below = int. Above x and below y = acc
    path_to_cov = path to directory with synteny
    
    return: 
    df row. Per column: mean sharedness of ref-acc with genome
    """
    #read df
    acc_df = read_df(f'{path_to_cov}/{genome}_covdf.csv')
    #get acc windows
    acc_df = acc_df[(acc_df['count'] >= above) & (acc_df['count'] <= below)]
    #per acc window determine mean number of windows shared per genome
    synt_regions = (acc_df.iloc[:,7:-1] >= 1).mean()
    return synt_regions

def count_syntenic(genome_list, path_to_cov, above, below):#, indir, ext):
    """count syntenic regions between two isolates
    A cutoff can specify regions present in at least x isolates 
    Interested in all? use above 0, below len(list).
    
    return:
    dataframe. Row = reference, col = comparison
    """
    syntenic_regions = []
    for ref in list(genome_list):
        try:
            #get acc regions and sharedness per genome for ref
            temp_sregions = \
            get_shared_per_genome(ref, above, below, path_to_cov)
            if len(syntenic_regions) != 0:
                #if synt region exists, concetanate temp region
                syntenic_regions = \
                    pd.concat((pd.DataFrame(syntenic_regions), \
                    pd.DataFrame(temp_sregions)), axis = 1)
                cols.append(ref)
            else:
                #if not yet syntenic regions, initiate
                syntenic_regions = temp_sregions
                cols = [ref]
        #if ref not present, print ref
        except KeyError:
            print(genome)
    syntenic_regions.columns = cols
    return syntenic_regions

if __name__ == '__main__':
    #read meta_data
    metadata = pd.read_csv('/home/anouk/SPP_FOCReSequencingMetadata.csv', sep = ',')
    #get genome_list
    genome_list = [x.split('_covdf')[0] for x in os.listdir(argv[1])]
    cutoff = len(genome_list) - 2
    #fill sizes dict
    genome_size_dict = fill_sizes_dict(genome_list, argv[1], cutoff)
    #plot sizes dict
    out_prefix = argv[3]
    plot_sizes(genome_size_dict, out_prefix)
    #get mean syntenic_dataframe, argv[1]: 'normal' argv[2]: sliding
    synt_reg = count_syntenic(genome_list, argv[1], 0, 0.8*len(genome_list))
    #plot syntenic dataframe as clustermap
    sns.clustermap(synt_reg.fillna(1), figsize= (40,40))
    #plt.show()
    plt.savefig(f'{out_prefix}_clustermap.pdf')
