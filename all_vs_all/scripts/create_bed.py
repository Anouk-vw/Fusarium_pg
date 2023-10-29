import pandas as pd
from sys import argv
import numpy as np

def read_covdf(covDf):
    """
    Get covDf
    """
    cov_pd = pd.read_csv(covDf)
    return cov_pd

def acc_regions(genome,cutoffs, path_to_covdf, category):
    #make_df
    cov_df = read_covdf(f'{path_to_covdf}{genome}_covdf.csv')
    #get accessory regions
    print('calculating regions')
    acc_regions = cov_df[(cov_df['count'] < cutoffs*0.8) & (cov_df['count'] > 1)][['Chrom', 'start', 'end']]
    """for i in set(list(bed_info['Chrom'])):
        tart = min(bed_info[bed_info['Chrom'] == i]['start'])
        end = max(bed_info[bed_info['Chrom']== i]['end'])
        if end - start > 100000:
            print(i,start,end)
    #bed_info.to_csv('accessory_bed/{}_accinfo.bed'.format(i), sep = '\t')"""
    return acc_regions

def get_consecutive(acc_regions):
    acc_regions_dict = {}
    #per chromsome with accessory regions, get longest stretch
    for i in set(list(acc_regions['Chrom'])):
        chrom = acc_regions[acc_regions['Chrom'] == i]
        #get difference between starts in chrom
        acc_starts = np.diff(chrom['start'])
        #add 5000 to make up for final element in list
        list(acc_starts).append(5000)
        consecutive = []
        temp = []
        for index, value in enumerate(acc_starts):
            #if value is below 50 000, the distance between acc regions is less than 50kb
            if value <= 50000:
                temp.append(index)
            #if index is the last item, append temp to consecutive
            if index == len(acc_starts) - 1:
                if len(temp) > 20:
                    consecutive.append(temp)
            #if value smaller than 50 000, break consecutive list and append temp to consecutive
            if value > 50000:
                #temp should at least contain 20 indexes = 5000*20 = 100kb consecutive acc
                if len(temp) > 20:
                    consecutive.append(temp)
                temp = []
        acc_regions_dict[i] = consecutive
    return acc_regions_dict

#for each genome get all consecutive acc regions and save them to a bed file (start stop)
def acc_bed(genome, cutoff, path_to_cov, category):
    bed_list = []
    ar = acc_regions(genome, cutoff, path_to_cov, category)
    ar_dict = get_consecutive(ar)
    for k,v in ar_dict.items():
        #if len(v) >= 1:
            #print(k,len(v))
        chrom_df = ar[ar['Chrom'] == k]
        for index_list in v:
            start_i = index_list[0]
            end_i = index_list[-1]
            start = chrom_df.iloc[start_i,1]
            end = chrom_df.iloc[end_i,2]
            bed_list.append([k,start,end])
    if len(bed_list) >0:
        bed_df = pd.DataFrame(bed_list)
        print(f'writing to {genome}_accinfo.bed')
        bed_df.to_csv('{}_accinfo.bed'.format(genome), sep = '\t', index = False, header = False)
    else:
        print(genome, bed_list)
    return

if __name__ == '__main__':
    print('hello')
    acc_bed(argv[1], 69, argv[2], 'accessory')
