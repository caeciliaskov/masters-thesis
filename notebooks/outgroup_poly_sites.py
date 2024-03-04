#A. Import
import sys
import vcf_processing_functions as vpf
import pandas as pd
import zarr
import numpy as np
import allel
from variables import *

#B. Variables
dataset = sys.argv[1]
datazar = sys.argv[2]
chrom   = sys.argv[3]
callset = zarr.open_group("{sandbox}/zar/{datazar}/{datazar}_chr{chrom}.zarr".format(sandbox = sandbox, datazar = datazar, chrom = chrom), mode='r')

#C. Code
#C.1. Get the index of samples that form the outgroup set
outgroup     = np.loadtxt("{results}/outgroup_{dataset}.txt".format(results = results, dataset = dataset), dtype=str)
samples      = callset["{chrom}/samples".format(chrom = chrom)][:]
outgroup_idx = np.array([i for i in range(len(samples)) if samples[i] in outgroup])

#C.2. Get the filters per dataset and the ancestral and derived alleles map
all_filt, max_min_anc_filt, anc_der_map, ref_alt_alleles, anc_der_map_correct = vpf.get_filters_and_map(chrom = chrom, dataset = dataset, callset = callset, ingroup_idx = outgroup_idx, vervose = True)

#C.3. Get the variable positions in the outgroup
count_alleles_outgroup = (allel.GenotypeDaskArray(callset['{chrom}/calldata/GT'.format(chrom = chrom)])
                               .take(outgroup_idx, axis = 1)     
                               .compress(all_filt, axis = 0)
                               .compress(max_min_anc_filt, axis = 0)
                               .map_alleles(anc_der_map)
                               .map_alleles(anc_der_map_correct)
                               .count_alleles()
                               .compute())

variant_loci_outgroup = count_alleles_outgroup.is_variant()


#C.4. Get the genomic coordinates of those variable positions in the outgroup and output them in a file
(pd.DataFrame({"dat" : dataset,
               "chr" : chrom,
               "pos" : allel.SortedIndex(callset["{chrom}/variants/POS".format(chrom = chrom)])[all_filt][max_min_anc_filt][variant_loci_outgroup],
               "anc" : ref_alt_alleles[max_min_anc_filt][np.argwhere(anc_der_map == 0)[:, 0], np.argwhere(anc_der_map == 0)[:, 1]][variant_loci_outgroup],
               "aco" : count_alleles_outgroup[variant_loci_outgroup, 0],
               "der" : ref_alt_alleles[max_min_anc_filt][np.argwhere(anc_der_map == 1)[:, 0], np.argwhere(anc_der_map == 1)[:, 1]][variant_loci_outgroup],
               "dco" : count_alleles_outgroup[variant_loci_outgroup, 1],
               "ref" : anc_der_map[variant_loci_outgroup, 0] == 0})
    .to_csv("{sandbox}/out/{dataset}/{dataset}_chr{chrom}_der.txt".format(sandbox = sandbox, dataset = dataset, chrom = chrom), 
                                                                          sep='\t', header=False, index=False))

#C.5. Get the polymorphic positions in the outgroup (not fixed for the anc or der)
count_alleles_outgroup = (allel.GenotypeDaskArray(callset['{chrom}/calldata/GT'.format(chrom = chrom)])
                               .take(outgroup_idx, axis = 1)     
                               .compress(all_filt, axis = 0)
                               .compress(max_min_anc_filt, axis = 0)
                               .map_alleles(anc_der_map)
                               .map_alleles(anc_der_map_correct)
                               .count_alleles()
                               .compute())

variant_loci_outgroup = count_alleles_outgroup.is_segregating()

#C.6. Get the genomic coordinates of those variable positions in the outgroup and output them in a file
(pd.DataFrame({"dat" : dataset,
               "chr" : chrom,
               "pos" : allel.SortedIndex(callset["{chrom}/variants/POS".format(chrom = chrom)])[all_filt][max_min_anc_filt][variant_loci_outgroup],
               "anc" : ref_alt_alleles[max_min_anc_filt][np.argwhere(anc_der_map == 0)[:, 0], np.argwhere(anc_der_map == 0)[:, 1]][variant_loci_outgroup],
               "aco" : count_alleles_outgroup[variant_loci_outgroup, 0],
               "der" : ref_alt_alleles[max_min_anc_filt][np.argwhere(anc_der_map == 1)[:, 0], np.argwhere(anc_der_map == 1)[:, 1]][variant_loci_outgroup],
               "dco" : count_alleles_outgroup[variant_loci_outgroup, 1],
               "ref" : anc_der_map[variant_loci_outgroup, 0] == 0})
    .to_csv("{sandbox}/out/{dataset}/{dataset}_chr{chrom}_seg.txt".format(sandbox = sandbox, dataset = dataset, chrom = chrom), 
                                                                          sep='\t', header=False, index=False))

# #C.7. Get number of SNPs (het = 1, hom der = 2) per individual per site
# g = (allel.GenotypeDaskArray(callset['{chrom}/calldata/GT'.format(chrom = chrom)])
#                               .take(outgroup_idx, axis = 1)     
#                               .compress(all_filt, axis = 0)
#                               .compress(max_min_anc_filt, axis = 0)
#                               .map_alleles(anc_der_map)
#                               .compute())

# count_SNPs_pedind = g.is_hom_alt()*2 + g.is_het()*1
    

# #C.8. Get the genomic coordinates of those positions, and sum the count of SNPs per individual in windows
# pos_filt = allel.SortedIndex(callset["{chrom}/variants/POS".format(chrom = chrom)])[all_filt][max_min_anc_filt]
# wei      = np.loadtxt("{sandbox}/wei/chr{chrom}.txt".format(sandbox = sandbox, chrom = chrom), usecols=[1, 2]) 

# count_SNPs_pedind_window = []

# for i in range(wei.shape[0]): 
#     count_SNPs_pedind_window.append(count_SNPs_pedind[(pos_filt > wei[i, 0]) * (pos_filt <= wei[i, 0]+window_len)].sum(axis = 0).tolist())
# count_SNPs_pedind_window = np.array(count_SNPs_pedind_window)
# #C.9. Write the output
# np.savetxt("{sandbox}/out/{dataset}/{dataset}_chr{chrom}_snp.txt".format(sandbox = sandbox, dataset = dataset, chrom = chrom),
#            count_SNPs_pedind_window, delimiter = "\t")

