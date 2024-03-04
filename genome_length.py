import pandas as pd
import zarr
import allel
from variables import *
import numpy as np
import vcf_processing_functions as vpf
import sys

dataset = 'HGDP'

outgroup = np.loadtxt("{project}/people/moi/results/outgroup_{dataset}.txt".format(project = project, dataset = dataset), dtype=str)

pos_length_list = []
chrom_list = []


for chrom in lchromos:
    chrom_list.append(chrom)
    callset = zarr.open_group("/home/clsj/GenerationInterval/people/moi/sandbox/zar/{dataset}/{dataset}_chr{chrom}.zarr".format(dataset = dataset,chrom = chrom), mode='r')
    samples = callset["{chrom}/samples".format(chrom = chrom)][:]
    ingroup_idx = np.array([i for i in range(len(samples)) if samples[i] not in outgroup])
    all_filt, max_min_anc_filt, anc_der_map, ref_alt_alleles, anc_der_map_correct = vpf.get_filters_and_map(chrom = chrom, dataset = dataset, callset = callset, ingroup_idx = ingroup_idx, vervose = True)
    pos = callset["{chrom}/variants/POS".format(chrom = chrom)][:][all_filt][max_min_anc_filt]
    idx = allel.SortedIndex(pos)
    pos_length_list.append(idx[-1]-idx[0]+1)

d = {'chr': chrom_list,'length':pos_length_list}
df = pd.DataFrame(data=d)

df.to_csv("~/GenerationInterval/people/clsj/script_results/tables/genome_length_Human.txt", index = False, header = True, sep = '\t')

