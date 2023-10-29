"""
python scripts/itterate_counts.py brocoli_2_noMINI/dir_step3/table_OGs_protein_counts.txt brocoli_2_noMINI/singletons_counts.txt itterate 
python ~/anouk2/pangenome/cactus_complete/genes/scripts/itterate_counts.py dir_step3/table_OGs_protein_counts.txt new_protein_singletons.txt itterate

"""
from sys import argv
import pandas as pd
import matplotlib.pyplot as plt
import random
import numpy as np

def get_singleton(single_list):
    single_dict = {}
    with open(single_list) as sl:
        for line in sl:
            genome, singletons = line.strip().split('\t')
            single_dict[genome] = singletons
    return single_dict

def calculate_size(dataframe, cutoff):
    COGs = dataframe[((dataframe != 0).sum(axis = 1) >= cutoff)].index.tolist()
    AOGs = dataframe[((dataframe != 0).sum(axis = 1) < cutoff)].index.tolist()
    #size of pangenome = total number of OGs present (=non zero)
    total = dataframe[((dataframe != 0).sum(axis = 1) != 0)].index.tolist()
    return len(COGs), len(AOGs), len(total)

def get_panSizes(genome_list, cTable, genome_names, sdict):
    """

    """
    c_OGs = []
    a_OGs = []
    u_OGs = []
    for i,genome in enumerate(genome_list):
        if i == 0:
            genomes = [genome_list[0]]
            cutoff = 1
            dataframe = cTable.iloc[:,genomes]
            core, acc, total = calculate_size(dataframe, cutoff)
            s_counts = int(sdict[genome_names[i]])
            c_OGs.append(core)
            a_OGs.append(total)
            u_OGs.append(total + s_counts)
        else:
            genomes = genome_list[0:i+1]
            cutoff = len(genomes)
            dataframe = cTable.iloc[:,genomes]
            s_counts = sum([int(sdict[genome_names[itt]]) for itt in range(0,i)])
            core, acc, total = calculate_size(dataframe, cutoff)
            c_OGs.append(core)
            a_OGs.append(total)
            u_OGs.append(total + s_counts)
    return c_OGs, a_OGs, u_OGs
    
def get_random(it, numGen):
    gene_orders = []
    for i in range(0,it):
        randomOrder = random.sample(range(0,numGen), numGen)
        gene_orders.append(randomOrder)
    return gene_orders

def multiple(orders, cTable, genome_names, sdict):
    c_list = []
    a_list = []
    u_list = []
    #print(orders)
    for i in orders:
        c, a, u = get_panSizes(i, cTable, genome_names, sdict)
        c_list.append(c)
        a_list.append(a)
        u_list.append(u)
    return c_list, a_list, u_list

def plot(core, u, a, output):
    #core
    means  =[np.mean([x[i] for x in core]) for i in range(len(core[0]))]
    errormax = [np.max([x[i] for x in core]) for i in range(len(core[0]))]
    errormin = [np.min([x[i] for x in core]) for i in range(len(core[0]))]
    errormax = np.subtract(means, errormax)
    errormin = np.subtract(errormin, means)
    yerr = np.transpose([[emax, emin] for emax,emin in zip(errormax, errormin)])
    # U + A = total
    meansu  =[np.mean([x[i] for x in u]) for i in range(len(u[0]))]
    errormaxu = [np.max([x[i] for x in u]) for i in range(len(u[0]))]
    errorminu = [np.min([x[i] for x in u]) for i in range(len(u[0]))]
    errormaxu = np.subtract(meansu,errormaxu)
    errorminu = np.subtract(errorminu, meansu)
    yerru = np.transpose([[emaxu, eminu] for emaxu,eminu in zip(errormaxu, errorminu)])
    #acc
    meansa  =[np.mean([x[i] for x in a]) for i in range(len(a[0]))]
    errormaxa = [np.max([x[i] for x in a]) for i in range(len(a[0]))]
    errormina = [np.min([x[i] for x in a]) for i in range(len(a[0]))]
    errormaxa = np.subtract(meansa,errormaxa)
    errormina = np.subtract(errormina, meansa)
    yerra = np.transpose([[emaxa, emina] for emaxa,emina in zip(errormaxa, errormina)])
    #plot
    plt.subplots(figsize=(15,10))
    #core
    plt.errorbar(x = range(len(core[0])), y = means, yerr= yerr, c = '#44803F', alpha = 0.7, fmt = '.')
    plt.scatter(x=list(range(len(means))), y = means, c= '#44803F', s = 60)
    #ototal
    plt.errorbar(x = range(len(u[0])), y = meansu, yerr= yerru, c = '#DBB679', alpha = 0.7, fmt ='.')
    plt.scatter(x=list(range(len(meansu))), y = meansu, c='#DBB679', s = 60)
    #Accessory only
    plt.errorbar(x = range(len(a[0])), y = meansa, yerr= yerra, c = '#FF5A33', alpha = 0.7, fmt ='.')
    plt.scatter(x=list(range(len(meansa))), y = meansa, c='#FF5A33', s = 60)
    plt.ylabel('Number of OG', fontsize = 42)
    plt.xlabel('Number of Genomes', fontsize = 42)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize = 24, rotation = 45)
    #plt.show()
    plt.grid(axis = 'y')
    plt.savefig(f'{output}.pdf', format ='pdf')
    return

def add_category(dataframe, cutoff):
    cat =[]
    for x in list((dataframe!= 0).sum(axis = 1)):
        x = int(x)
        if int(x) >= cutoff:
            cat.append('core')
        if x < cutoff and x >= cutoff-10:
            cat.append('softcore')
        if x < cutoff-10 and x>=2:
            cat.append('accessory')
        if x < 2:
            cat.append('unique')
    dataframe['category'] = cat
    return dataframe

def plot_CA_genome(catdf, sdict, names):
    core = list((catdf[catdf['category'] == 'core'] != 0).sum())[0:-1]
    softcore = list((catdf[catdf['category'] == 'softcore'] != 0).sum())[0:-1]
    accessory= list((catdf[catdf['category'] == 'accessory'] != 0).sum())[0:-1]
    unique = [int(sdict[x]) for x in names] #list(catdf[catdf['category'] == 'unique'].sum())[0:-1] + [int(sdict[x]) for x in names]
    print(unique)
    #print([(x,y) for x,y in zip(catdf.columns[:-1],names)])
    plt.subplots(figsize = (10,10))
    plt.bar(names, core)
    plt.bar(names, softcore, bottom = core)
    plt.bar(names, accessory, bottom = np.array(core) + np.array(softcore))
    plt.bar(names, unique, bottom = np.array(core) + np.array(softcore) + np.array(accessory))
    plt.xticks(rotation = 90)
    plt.show()

if __name__ == '__main__':
    extension = argv[3]
    s_genes = get_singleton(argv[2])
    tableCounts = pd.read_csv(argv[1], sep = '\t', index_col=0)
    tableCounts.columns = [x.replace(extension, '') for x in tableCounts.columns]
    random_itterations = get_random(10,69)
    gene_list = tableCounts.columns #[x.replace(extension, '') for x in tableCounts.columns]
    if 'itterate' in argv[4]:
        core_list, acc_list, u_list= multiple(random_itterations, tableCounts, gene_list, s_genes)
        plot(core_list, u_list, acc_list, argv[5])
    if 'size' in argv[4]:
        categories = add_category(tableCounts, 69)
        plot_CA_genome(categories, s_genes, gene_list)
