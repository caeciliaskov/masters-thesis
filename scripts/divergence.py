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

pos = callset["{chrom}/variants/POS".format(chrom = chrom)][:]

#########################################################
#     Sample individuals from selected populations      #
#########################################################

# Here I take the populations that I want and put them in a list (pop_list) and I create a new list to put the indices in for each population (pop_idx).
# For each population in pop_list, I then sample all the individuals from that population and take the indices and add them to the pop_idx list. Pop_idx is then a list of lists.
# I then flatten pop_idx to create just a list and then I sort it to use it in the dask array.

pop_list = ['Burusho','Hazara','Uygur','Bougainville','PapuanHighlands','PapuanSepik','Colombian','Karitiana','Surui','Pima', 'Maya','Yakut','Oroqen','Daur','Hezhen','Mongolian','Xibo','Japanese','NorthernHan','Naxi','Yi','Tu','Tujia','Cambodian','Dai','Lahu','Miao','Han','She','Adygei','Basque','Sardinian','BergamoItalian','Tuscan','French','Orcadian','Bedouin','Palestinian','Druze','Russian','Balochi','Brahui','Makrani','Kalash','Pathan','Sindhi','Mozabite']
pop_idx = []

for x in pop_list:
    data = metadata[metadata['pop'] == x]
    samples = data.sample(n=len(data), random_state=10)
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
#                    Create filters                     #
#########################################################

# Get the filters per dataset and the ancestral and derived alleles map

all_filt, max_min_anc_filt, anc_der_map, ref_alt_alleles, anc_der_map_correct = vpf.get_filters_and_map(chrom = chrom, dataset = dataset, callset = callset, ingroup_idx = idx_list, vervose = True)

#########################################################
#                     Count alleles                     #
#########################################################

# Here I create a dictonary for the allele counts. For each list number x in pop_idx, I create a key with the name of the population and then I create allele counts arrays.
    # With take() I only take the indices from the idx_list.
    # With compress() I apply the allele filter that filters for positions that are not included in the bed file ranges (acc_filt), PASS flag filter FILTER_ExcHet (pas_filt), indel positions is_snp (snp_filt), positions that have only first and second maximum frequency alleles (max_min_filt), excess heterozygosity (exc_filt), and HWE test (hwe_filt)
    # With compress() I filter for ancestry (max_min_anc_filt).

ac_list = {}
for x in range(len(pop_idx)):
    ac_list[pop_list[x]] = (allel.GenotypeDaskArray(callset["{chrom}/calldata/GT".format(chrom = chrom)])
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

i_list = []
j_list = []
d_list = []
for i in ac_list:
    for j in ac_list:
        i_list.append(i)
        j_list.append(j)
        d_list.append(allel.sequence_divergence(pos, ac_list[i], ac_list[j]))

#########################################################
#                       Make plot                       #
#########################################################

# The Dxy values are then saved in a dataframe together with the two populations compared for each value. The D values are then normalized using z-score, and the diagonal is manually changed to 0.
# The heat map is then made by manually setting the order for column and row.

# Save dataframe with Dxy values:
d = {'i': i_list, 'j': j_list, 'D': d_list}
df = pd.DataFrame(data=d)

# Normalize data using z-score:
df_sklearn = df.copy()
df_sklearn['D_norm'] = zscore(df_sklearn['D'])

# Change values in the diagonal to zero:
df_sklearn.loc[df["i"] == df['j'], "D_norm"] = 0

# Make plot:
result = df_sklearn.pivot(index='i',columns='j',values='D_norm')

column_order = ['Burusho','Hazara','Uygur','Bougainville','PapuanHighlands','PapuanSepik','Colombian','Karitiana','Surui','Pima', 'Maya','Yakut','Oroqen','Daur','Hezhen','Mongolian','Xibo','Japanese','NorthernHan','Naxi','Yi','Tu','Tujia','Cambodian','Dai','Lahu','Miao','Han','She','Adygei','Basque','Sardinian','BergamoItalian','Tuscan','French','Orcadian','Bedouin','Palestinian','Druze','Russian','Balochi','Brahui','Makrani','Kalash','Pathan','Sindhi','Mozabite']
row_order = ['Burusho','Hazara','Uygur','Bougainville','PapuanHighlands','PapuanSepik','Colombian','Karitiana','Surui','Pima', 'Maya','Yakut','Oroqen','Daur','Hezhen','Mongolian','Xibo','Japanese','NorthernHan','Naxi','Yi','Tu','Tujia','Cambodian','Dai','Lahu','Miao','Han','She','Adygei','Basque','Sardinian','BergamoItalian','Tuscan','French','Orcadian','Bedouin','Palestinian','Druze','Russian','Balochi','Brahui','Makrani','Kalash','Pathan','Sindhi','Mozabite']

table3 = result[column_order][row_order]
table4 = table3.reindex(row_order)

fig, ax = plt.subplots(figsize=(8, 6))
heat_map = sns.heatmap(table4, linewidth = 0, xticklabels=True,yticklabels=True, cmap = sns.color_palette("blend:#1B0E0A,#8A1418,#E41D26,#ED5623,#F8B40D,#F8EC3A,#FCF7C1", as_cmap=True))
plt.savefig('~/GenerationInterval/people/clsj/script_results/chr{chrom}.png'.format(chrom = chrom))





