import allel
import pandas as pd
import zarr
import numpy as np
import vcf_processing_functions as vpf
from variables import *
import sys
from scipy.stats import zscore
import matplotlib.pyplot as plt
from sklearn import preprocessing
import seaborn as sns
import os
import itertools

sys.path.insert(1, '/home/clsj/GenerationInterval/people/moi/scripts')

#########################################################
#                    Load variables                     #
#########################################################

data1 = sys.argv[1]
anc1 = sys.argv[2]
data2 = sys.argv[3]
anc2 = sys.argv[4]
mean_prob = sys.argv[5]
chrom = sys.argv[6]

#########################################################
#                      Load data                        #
#########################################################

def get_int_total_seq_length(pop1, pop2):
    return (pd.read_csv("{tables}/neanderthal_summary_table.txt".format(tables=tables),sep='\t',header=0)
                   .query('pop1 == "{}" & pop2 == "{}"'.format(pop1, pop2))["int_total_seq_length"]
                   .to_numpy()[0])

callset = zarr.open_group("{project}/people/moi/sandbox/zar/{dataset}/{dataset}_chr{chrom}.zarr".format(dataset = data1,chrom = chrom, project = project), mode='r')

samples = callset["{chrom}/samples".format(chrom = chrom)][:]

metadata = pd.read_table("~/GenerationInterval/people/moi/files/metadata.txt")
metadata = metadata[metadata['dat'] == data1]

samples_dataset = pd.DataFrame({'ind': samples})
metadata = pd.merge(samples_dataset, metadata, how='inner', on = 'ind')

#########################################################
#                     Create filter                     #
#########################################################

outgroup = np.loadtxt("{project}/people/moi/results/outgroup_{dataset}.txt".format(project = project, dataset = data1), dtype=str)
ingroup_idx = np.array([i for i in range(len(samples)) if samples[i] not in outgroup])

all_filt, max_min_anc_filt, anc_der_map, ref_alt_alleles, anc_der_map_correct = vpf.get_filters_and_map(chrom = chrom, dataset = data1, callset = callset, ingroup_idx = ingroup_idx, vervose = True)

pos = callset["{chrom}/variants/POS".format(chrom = chrom)][:][all_filt][max_min_anc_filt]


#########################################################
#                  Make genotype array                  #
#########################################################

g = (allel.GenotypeDaskArray(callset["{chrom}/calldata/GT".format(chrom = chrom)])
                                    .take(samples_dataset.index.to_list(),axis=1)
                                    .compress(all_filt, axis = 0)
                                    .compress(max_min_anc_filt, axis = 0)
                                    .map_alleles(anc_der_map)
                                    .map_alleles(anc_der_map_correct)
                                    .compute())

#########################################################
#                      Load data                        #
#########################################################

for pop1, pop2 in itertools.combinations(lpopulat, 2):
    data = pd.read_csv("{}/int_{}_{}_{}_{}_{}_{}_{}.txt".format(intersection,data1,pop1,anc1,data2,pop2,anc2,mean_prob), sep='\t', header=None)
    data.columns = ['chr', 'start', 'end', 'ind1', 'ind2']

#########################################################
#                     Get positions                     #
#########################################################

    ranges = []
    for c in data.values:
        if c[0] == 'chr{chrom}'.format(chrom=chrom):
            ranges.append([c[1],c[2]])

    idx = allel.SortedIndex(pos)

    ranges = np.array(ranges)

    starts = ranges[:, 0]
    stops = ranges[:, 1]

    frag_pos = []

    for i in range(len(ranges)):
        try:
            loc = idx.locate_range(starts[i], stops[i])
            pos_in = idx[loc]
            frag_pos.append(list(pos_in))
        except KeyError:
            print('There was a KeyError at {} {}'.format(starts[i], stops[i]))
            frag_pos.append([])
        except Exception as e: 
            print(e)

    pos_idx_flat = (pd.DataFrame({"idx" : np.arange(idx.shape[0]), 
                              "pos" : idx})
                      .merge(pd.DataFrame({"pos" : np.concatenate(frag_pos)})))["idx"].to_numpy()

    pos_idx = []
    prev_i = 0
    for x in frag_pos:
        i = len(x)
        pos_idx.append(pos_idx_flat[prev_i:i])
        prev_i = i


#########################################################
#               Make allele counts array                #
#########################################################

    ac_list = {}
    pop_list = [pop1, pop2]
    for i in range(len(pop_list)):
        fragment_procedence_sorted = data.query('chr == "chr{}"'.format(chrom))["ind{}".format(i+1)].to_numpy()
        number_of_snps_in_each_frag = np.array([len(pos_idx[i]) for i in range(len(pos_idx))])
        variant_procedence_sorted = np.repeat(fragment_procedence_sorted, number_of_snps_in_each_frag)
        variant_procedence_index_sorted = np.array([np.where(samples == variant_procedence_sorted[i])[0][0] for i in range(variant_procedence_sorted.shape[0])])
        artificial_variant_genome_p0 = g.take(pos_idx_flat, axis=0)[np.arange(variant_procedence_index_sorted.shape[0]), variant_procedence_index_sorted]
        artificial_variant_genome_p0 = artificial_variant_genome_p0.reshape(artificial_variant_genome_p0.shape[0], 1, artificial_variant_genome_p0.shape[1])
        ac = allel.GenotypeArray(artificial_variant_genome_p0).count_alleles()
        ac_list[pop_list[i]] = ac

#########################################################
#                     Calculate Dxy                     #
#########################################################



    i_list = []
    j_list = []
    d_list = []
    mpd_list = []
    n_bases = []
    for i in ac_list:
        for j in ac_list:
            i_list.append(i)
            j_list.append(j)
            d_list.append(allel.sequence_divergence(np.arange(ac_list[i].shape[0]), ac_list[i], ac_list[j]))
            mpd = allel.mean_pairwise_difference_between(ac_list[i], ac_list[j])
            mpd_sum = np.sum(mpd)
            if i == j:
                n_bases.append('0')
            else:
                n_bases.append(get_int_total_seq_length(pop1, pop2))
            mpd_list.append(mpd_sum)


#########################################################
#                      Save to file                     #
#########################################################

    d = {'i': i_list, 'j': j_list, 'D': d_list,'mpd':mpd_list,'n_bases':n_bases}
    df = pd.DataFrame(data=d)
    df.to_csv("{}/div_{}_{}_{}_{}_{}_{}_{}_chr{}.txt".format(divergence, data1, pop1, anc1, data2, pop2, anc2, mean_prob, chrom), sep='\t', index=False, header=None)

    with open(os.path.join(completed, 'div_{}_{}_{}_{}_{}_{}_{}_chr{}.DONE'.format(data1, pop1, anc1, data2, pop2, anc2, mean_prob, chrom)), 'w') as fp:
        pass