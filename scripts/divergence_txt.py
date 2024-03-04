import allel
import pandas as pd
import seaborn as sns
import random
import zarr
import numpy as np
import vcf_processing_functions as vpf
from variables import *
import sys
import matplotlib.pyplot as plt
from sklearn import preprocessing
from scipy.stats import zscore
import scipy

sys.path.insert(1, '/home/clsj/GenerationInterval/people/moi/scripts')

#########################################################
#                      Load data                        #
#########################################################

dataset = sys.argv[1]
chrom = sys.argv[2]
callset = zarr.open_group("/home/clsj/GenerationInterval/people/moi/sandbox/zar/{dataset}/{dataset}_chr{chrom}.zarr".format(dataset = dataset,chrom = chrom), mode='r')

samples = callset["{chrom}/samples".format(chrom = chrom)][:]

metadata = pd.read_table("~/GenerationInterval/people/moi/files/metadata.txt")
metadata = metadata[metadata['dat'] == dataset]

samples_dataset = pd.DataFrame({'ind': samples})
metadata = pd.merge(samples_dataset, metadata, how='inner', on = 'ind')

#########################################################
#                    Create filters                     #
#########################################################

# Get the filters per dataset and the ancestral and derived alleles map

outgroup = np.loadtxt("{project}/people/moi/results/outgroup_{dataset}.txt".format(project = project, dataset = dataset), dtype=str)
ingroup_idx = np.array([i for i in range(len(samples)) if samples[i] not in outgroup])

all_filt, max_min_anc_filt, anc_der_map, ref_alt_alleles, anc_der_map_correct = vpf.get_filters_and_map(chrom = chrom, dataset = dataset, callset = callset, ingroup_idx = ingroup_idx, vervose = True)

pos = callset["{chrom}/variants/POS".format(chrom = chrom)][:][all_filt][max_min_anc_filt]

#########################################################
#     Sample individuals from selected populations      #
#########################################################

# Here I take the populations that I want and put them in a list (pop_list) and I create a new list to put the indices in for each population (pop_idx).
# For each population in pop_list, I then sample all the individuals from that population and take the indices and add them to the pop_idx list. Pop_idx is then a list of lists.
# I then flatten pop_idx to create just a list and then I sort it to use it in the dask array.

pop_idx = []

for c in lpopulat:
    data = metadata[metadata['pop'] == c]
    samples = data.sample(n=len(data))
    idx = samples.index.tolist()
    pop_idx.append(idx)

idx_list = [item for sublist in pop_idx for item in sublist]
idx_list.sort()

#########################################################
#        Create list with subpopulation indices         #
#########################################################

# For each population in pop_idx, I create a new list. This list l will consist of the new index in the idx_list for each individual index in pop_idx. Each new list l with the new indices for each population, will then be added to a list called subpops.

subpops = []

for i in pop_idx:
    l = []
    for j in i:
        l.append(idx_list.index(j))
    subpops.append(l)    


#########################################################
#                     Count alleles                     #
#########################################################

# Here I create a dictonary for the allele counts. For each list number x in pop_idx, I create a key with the name of the population and then I create allele counts arrays.
    # With take() I only take the indices from the idx_list.
    # With compress() I apply the allele filter that filters for positions that are not included in the bed file ranges (acc_filt), PASS flag filter FILTER_ExcHet (pas_filt), indel positions is_snp (snp_filt), positions that have only first and second maximum frequency alleles (max_min_filt), excess heterozygosity (exc_filt), and HWE test (hwe_filt)
    # With compress() I filter for ancestry (max_min_anc_filt).

ac_list = {}
for x in range(len(pop_idx)):
    ac_list[lpopulat[x]] = (allel.GenotypeDaskArray(callset["{chrom}/calldata/GT".format(chrom = chrom)])
                                .take(idx_list, axis=1)
                                .compress(all_filt, axis = 0)
                                .compress(max_min_anc_filt, axis = 0)
                                .map_alleles(anc_der_map)
                                .map_alleles(anc_der_map_correct)
                                .count_alleles(subpop = np.array(subpops[x]))
                                .compute())


#########################################################
#                     Calculate Dxy                     #
#########################################################

# Calculating pairwise Dxy using sequence_divergence(). For each allele counts array in ac_list, the array will be compared pairwise to all other arrays and a Dxy value will be calculated and added to d_list.

pop1_list = []
pop2_list = []
mpd_list = []
for a in ac_list:
    for b in ac_list:
        pop1_list.append(a)
        pop2_list.append(b)
        mpd = allel.mean_pairwise_difference_between(ac_list[a], ac_list[b], fill=0)
        mpd_sum = np.sum(mpd)
        mpd_list.append(mpd_sum)


#########################################################
#                     Save dataframe                    #
#########################################################

# The Dxy values are then saved in a dataframe together with the two populations compared for each value. The D values are then normalized using z-score, and the diagonal is manually changed to 0.
# The heat map is then made by manually setting the order for column and row.

# Save dataframe with Dxy values:
d = {'pop1': pop1_list, 'pop2': pop2_list, 'mpd_human': mpd_list}
df = pd.DataFrame(data=d)

df.to_csv("~/GenerationInterval/people/clsj/script_results/tables/div_Human_chr{chrom}.txt".format(chrom = chrom), index = False, header = False, sep = '\t')

