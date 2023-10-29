def init_df(indir, genome_id,covcutoff, ext, sep):
    """
    Read all windows files as output by nucmer snakemake
    input: 
    genome = name of 'reference' genome
    ext = ending of the files to read ('sliding.bed')
    
    output:
    dataframe with windows of reference, coverage per genome per window
    """
    cov_dfs = []
    headers = []
    inputFiles = [x for x in os.listdir(indir) if x.startswith('{}{}'.format(genome_id, sep)) and x.endswith(ext)]
    for f in inputFiles:
        #print(f)
        df = pd.read_csv('{}{}'.format(indir, f), sep = '\t', header = None)
        if len(cov_dfs) > 0:
            cov_dfs.append(df[6])
            headers.append(f.split(sep)[1].split(ext)[0])
        else:
            cov_dfs.append(df)
            headers.extend(['Chrom', 'start', 'end', 'counts', 'bp', 'total', f.split(sep)[1].split(ext)[0]])
    cov_df = pd.concat(cov_dfs, axis = 1)
    cov_df = pd.concat((cov_df, (cov_df.iloc[:,6:,] >= covcutoff).sum(axis = 1)), axis = 1) #(cov_df.iloc[:,6:-1].mean(axis = 1))), axis = 1)
    headers.append('count')
    cov_df.columns = headers
    return cov_df
    
def save_df(indir, genome, ext, sep, cutoff, output):
    """
    Save dataframe
    """
    print(indir, genome , cutoff, ext, sep, output)
    df = init_df(indir, genome , cutoff, ext, sep)
    df.to_csv(output)
