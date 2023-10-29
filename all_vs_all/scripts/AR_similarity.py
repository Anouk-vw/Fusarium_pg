"""
Script to calculate pairwise similarity
cov_df > AR similarity
"""

def create_df(genome_id,covcutoff,ext,indir):
    cov_dfs = []
    headers = []
    inputFiles = [x for x in os.listdir(indir) if x.startswith('{}windows'.format(genome_id)) and x.endswith(ext) and not 'C068' in x]
    for f in inputFiles:
        df = pd.read_csv('{}{}'.format(indir, f), sep = '\t', header = None)
        if len(cov_dfs) > 0:
            cov_dfs.append(df[6])
            headers.append(f.split('windows')[1].split(ext)[0])
        else:
            cov_dfs.append(df)
            headers.extend(['Chrom', 'start', 'end', 'counts', 'bp', 'total', f.split('windows')[1].split(ext)[0]])
    cov_df = pd.concat(cov_dfs, axis = 1)
    cov_df.columns = headers
    cov_df = pd.concat((cov_df, (cov_df.iloc[:,6:,] >= covcutoff).sum(axis = 1)), axis = 1) #(cov_df.iloc[:,6:-1].mean(axis = 1))), axis = 1)
    headers.append('count')
    cov_df = pd.concat((cov_df, (cov_df.iloc[:,6:-1]).mean(axis = 1)), axis = 1) #(cov_df.iloc[:,6:-1].mean(axis = 1))), axis = 1)
    headers.append('synt ratio')
    cov_df.columns = headers
    return cov_df

def fill_df(genomes):
	"""
	genomes = list of genomes to include
	"""
	cutoff= len(genomes) - 2
	genomes_dict = {}
	for g in list(genomes):
		genomevs_df = create_df(g,0.8, ext, indir)
		genomes_dict = fill_core_dict(genomevs_df, genomes_dict, g, cutoff)
	return genomes_dict

def get_acc(genome, above, below, indir, ext):
    acc_df = create_df(genome, 0.8, ext, indir)
    acc_df = acc_df[(acc_df['count'] < below) & (acc_df['count'] > above)]
    return acc_df

def get_shared_per_genome(genome, above, below, indir, ext):
    """
    Return mean of shared regions between two isolates
    """
    #read df
    acc_df = create_df(genome, 0.8, ext, indir)
    #select regions based on sharedness
    acc_df = acc_df[(acc_df['count'] < below) & (acc_df['count'] > above)]
    #For these regions get mean sharedness per isolate
    synt_regions = (acc_df.iloc[:,6:-2] >= 0.8).mean()
    return synt_regions

def count_syntenic(genome_list, above, below, indir, ext):
    """count syntenic regions between two isolates
    A cutoff can be used to specify only looking into regions present in at least x isolates
    Almost unqiue? use below 3 and above 0. 
    All? use above 0, below len(list).
    Core? use above 80%
    Accessory? From 0 - 80%
    """
    #cols = []
    syntenic_regions = []
    #loop of genomes
    for g in list(genome_list):
        try:
            #get acc regions and sharedness per genome for g
            temp_sregions = get_shared_per_genome(g, above, below, indir, ext)
            if len(syntenic_regions) != 0:
                #concetanate temp region (if df is initialized)
                syntenic_regions = pd.concat((pd.DataFrame(syntenic_regions), pd.DataFrame(temp_sregions)), axis = 1)
                cols.append(g)
            else:
                #if not yet syntenic regions, initiate
                syntenic_regions = temp_sregions
                cols = [g]
        #if g not present, print g
        except KeyError:
            print(genome)
    syntenic_regions.columns = cols
    return syntenic_regions
    
def plot_ARperc(meta, synt_reg):
    """
    Plot boxplot with pairwise values between group of isoaltes
    """
    toplot_species = []
    toplot_labels = ['Tr4.R1', 'Host', 'Species', 'genotype', 'clade', 'location']
    for cat in toplot_labels:
        print(cat)
        temp= []
        for option in set(meta[cat]):
            cat_ids = list(meta[meta[cat] == option].IsolateCode)
            cat_ids = [x for x in cat_ids if x in synt_reg.index]
            if len(synt_reg.loc[cat_ids, cat_ids]) > 0:
                temp.extend(synt_reg.loc[cat_ids, cat_ids].stack().values)
        print(np.median(temp))
        toplot_species.append(temp)

    plt.violinplot(toplot_species, showmedians = True)
    labels = ['Race', 'Banana Host', 'Species', 'VCG', 'clade', 'Location']
    plt.xticks(range(1,len(toplot_labels)+1), labels, rotation = 90)
    return

    
if __name__ == '__main__':
    indir = 'nuc_synt/'
    ext = '_synt_5000.bed'
    genomes = set([x.split('windows')[0] for x in os.listdir(indir) if x.endswith(ext) and 'windows' in x and not 'C068' in x])
    
    #determine percentage of shared regions
    synt_regions = count_syntenic(genomes, 0,len(genomes)*0.8, indir, '_synt_5000.bed')
    
    #per group of isolates, plot all pairwise similarities as boxplot
    plot_ARperc(meta, synt_regions)
    
    #provide phylogenetic ordering of strains
    ordered_species_fob = ['36102','BRIP62280','BRIP63632','Leb1-2C',\
    'BRIP65062','M3','M2','M4','M5','S1B8','Race4','Pak1.1A', 'TR4_II5',\
    'II_5','M1','Eden','ISR1','ISR5','JV11','JV14','Laos','Vietnam',\
    'Myanmar','Col17','Col2','Col4','Indo8','Phi2-6C','9_KT06-A1',\
    'NRRL_36107','NRRL_36103','Bif_04','CR1.1','FocST4-98','FocP1',\
    'NRRL_36101','NRRL_36112','NRRL_36110','C058','C187','C082','C176',\
    'C177','C081','C192','C135','C090','NRRL_36117','NRRL_36108',\
    'P41b','NRRL_36118','19_KB07-B','56_JB09-A2','NRRL_36120',\
    'NRRL_36115','Mal43','NRRL_36116','F9129','15_KTG06-B2','Race1',\
    'Foc8','Foc16','P26c','P20a','Phi6.6a','Indo110']
    
    #order calculated similarity
    synt_regions_ordered = synt_regions.loc[ordered_species, ordered_species]
    
    #plot similarity heatmap
    sns.heatmap(synt_regions_ordered.fillna(1))
    plt.savefig('fsp_sharedAGRs_heatmap.pdf')
