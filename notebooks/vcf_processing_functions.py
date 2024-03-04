import allel
import numpy as np
from variables import *


def get_accessibility_mask_per_chrom(chrom):
    with open(accessmask, "r") as file:
        starts = []
        stops  = []
        for line in file:
            c, s, e = line.strip().split()
            if c == "chr{}".format(chrom):
                starts.append(int(s))
                stops.append(int(e))
    return np.array(starts), np.array(stops)

def read_ancestral_fasta(chrom, pos):
    j = 0
    i = 0
    anc = []
    with open("{data}/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_{chrom}.fa".format(data = data, chrom = chrom), "r") as file:
        for z, line in enumerate(file):
            if z:
                l = line.strip().upper()
                #j += len(l)
                if i >= len(pos):
                    break
                while pos[i] < j+len(l):
                    anc.append(l[pos[i]-1-j])
                    i += 1
                    if i >= len(pos):
                        break
                j += len(l)
                    
    return np.array(anc)

def maxmin_alleles(chrom, callset):
    #Total allele counts genotyped
    tot_alleles_count  = callset["{}/variants/AN".format(chrom)][:]
    #Allele counts for alternative alleles
    alt_alleles_count  = callset["{}/variants/AC".format(chrom)][:]
    alt_alleles_count[alt_alleles_count == -1] = 0
    #Allele counts for reference alleles
    ref_alleles_count  = tot_alleles_count-np.sum(alt_alleles_count, axis = 1)
    #Allele counts for all alleles together in a single array
    alleles_count      = np.concatenate((ref_alleles_count.reshape(-1, 1), alt_alleles_count), axis = 1)
    alleles_count_copy = np.copy(alleles_count)
    #Which is the allele with the maximum frequency?
    max_idx = np.argmax(alleles_count_copy, axis = 1)
    #Which is the allele with the second maximum frequency?
    alleles_count_copy[np.arange(alleles_count_copy.shape[0]), max_idx] = 0
    min_idx = np.argmax(alleles_count_copy, axis = 1)
    #max and min alleles index in a single array
    max_min_idx = np.concatenate((max_idx.reshape(-1, 1), min_idx.reshape(-1, 1)), axis = 1)
    #Positions that have only a First and second maximum frequency alleles
    alleles_count_copy[np.arange(alleles_count_copy.shape[0]), min_idx] = 0
    max_min_filt = np.prod(alleles_count[np.arange(alleles_count.shape[0]), min_idx].reshape(-1, 1) > alleles_count_copy, axis = 1).astype(np.bool_)
    return max_idx, min_idx, max_min_idx, max_min_filt

def pop_specific_filt(chrom, callset, stat, pops, max_idx, min_idx, threshold = 1e-6):
    filt = np.full(max_idx.shape, True)
    for pop in pops:
        #Get the array with the stat that we want to check specific for each population
        arr_orig = callset["{}/variants/{}_{}".format(chrom, stat, pop)][:]
        for maxmin in [max_idx, min_idx]:
            #Since the HWE and ExcHet are not performed for the reference allele, we want to check only positions in which the Max or min allele are not the reference allele
            maxmin_idx_no_ref      = maxmin != 0
            #From the original array with the test stat, we take take those positions that correspond to those alleles that are not the reference allele and check that the test is < than a certain value
            stat_filt                   = arr_orig[maxmin_idx_no_ref][np.arange(maxmin_idx_no_ref.sum()), maxmin[maxmin_idx_no_ref]-1] > threshold 
            #We keep this information for every population
            filt[maxmin_idx_no_ref] = (filt[maxmin_idx_no_ref])*stat_filt
    return filt

def get_anc_der_map(chrom, callset, all_filt, max_min_idx):
    max_min_idx      = max_min_idx[all_filt]
    pos              = allel.SortedIndex(callset["{}/variants/POS".format(chrom)])[all_filt]
    anc_pos          = read_ancestral_fasta(chrom, pos)
    anc_filt         = (anc_pos == "A")+(anc_pos == "C")+(anc_pos == "G")+(anc_pos == "T")
    ref_alt_alleles  = np.concatenate(((callset["{}/variants/REF".format(chrom)][:][all_filt]).reshape(-1, 1), callset["{}/variants/ALT".format(chrom)][:][all_filt]), axis = 1)
    max_min_alleles  = ref_alt_alleles[np.repeat(np.arange(ref_alt_alleles.shape[0]), max_min_idx.shape[1]), max_min_idx.reshape(-1)].reshape(-1, 2)
    max_min_anc_bool = max_min_alleles == np.repeat(anc_pos, 2).reshape(-1, 2)
    max_min_anc_filt = max_min_anc_bool.sum(axis = 1).astype("bool_")

    anc_idx = max_min_idx[max_min_anc_filt][max_min_anc_bool[max_min_anc_filt]]
    der_idx = max_min_idx[max_min_anc_filt][~max_min_anc_bool[max_min_anc_filt]]

    anc_der_map = np.full((max_min_anc_filt.sum(), ref_alt_alleles.shape[1]), -1)
    anc_der_map[np.arange(anc_der_map.shape[0]), anc_idx] = 0
    anc_der_map[np.arange(anc_der_map.shape[0]), der_idx] = 1
    
    return anc_filt, max_min_anc_filt, anc_der_map, ref_alt_alleles

def get_anc_der_map_correct(chrom, callset, ingroup_idx, all_filt, max_min_anc_filt, anc_der_map):
    count_alleles = (allel.GenotypeDaskArray(callset['{chrom}/calldata/GT'.format(chrom = chrom)])
                                    .take(ingroup_idx, axis = 1)     
                                    .compress(all_filt, axis = 0)
                                    .compress(max_min_anc_filt, axis = 0)
                                    .map_alleles(anc_der_map)
                                    .count_alleles()
                                    .compute())

    print("If there is a warning in the next line, it is because there are positions that there is no individuals called, so when the allele frequency is calculated, there is a division by 0")
    high_freq  = count_alleles[:, 1]/count_alleles.sum(axis = 1) > 0.95
    anc_is_alt = anc_der_map[:, 0] != 0


    anc_der_map_correct = np.tile(np.array([0, 1, -1, -1, -1, -1, -1]), anc_der_map.shape[0]).reshape(anc_der_map.shape)

    anc_der_map_correct[high_freq*anc_is_alt, 0] = 1
    anc_der_map_correct[high_freq*anc_is_alt, 1] = 0

    return anc_der_map_correct

def get_filters_and_map(chrom, dataset, callset, ingroup_idx, vervose = True):
    #A. BUILD FILTERS
    #A.1. Filter by Accessibility mask
    #A.1.1. Get accessibility bed file start and stop coordinates
    starts, stops = get_accessibility_mask_per_chrom(chrom)
    #A.1.2. Filter positions that are not included in the bed file ranges (bool)
    acc_filt = allel.SortedIndex(callset["{}/variants/POS".format(chrom)]).locate_ranges(starts, stops, strict = False)
    #A.2. PASS flag filter (bool)
    if dataset in ["1KGP", "ABOR", "VANU", "SGDP", "IGDP"]:
        pas_filt = callset["{}/variants/FILTER_PASS".format(chrom)][:]
    if dataset == "HGDP":
        pas_filt = ~callset["{}/variants/FILTER_ExcHet".format(chrom)][:]
    if dataset == "AYTA":
        pas_filt = ~callset["{}/variants/FILTER_LowQual".format(chrom)][:]
    #A.3. Filter out indel positions (bool)
    if dataset in ["SGDP"]:
        snp_filt = callset["{}/variants/numalt".format(chrom)][:] == 1
    else:
        snp_filt = callset["{}/variants/is_snp".format(chrom)][:]
    #A.4. Check the maximum frequency and second maximum frequency alleles for all positions
    if dataset in ["ABOR", "VANU", "SGDP", "AYTA", "IGDP"]:
        max_idx      = np.full(acc_filt.shape, 0)
        min_idx      = np.full(acc_filt.shape, 1)
        max_min_idx  = np.concatenate((max_idx.reshape(-1, 1), min_idx.reshape(-1, 1)), axis = 1)
        max_min_filt = np.full(acc_filt.shape, True)
    else:
        max_idx, min_idx, max_min_idx, max_min_filt = maxmin_alleles(chrom, callset)
    #A.5. Filter for excess of Heterozygosity and HWE test
    if dataset == "1KGP":
        exc_filt = pop_specific_filt(chrom, callset, "ExcHet", ["EUR", "AFR", "AMR", "SAS", "EAS"], max_idx, min_idx, 1e-6)
        hwe_filt = pop_specific_filt(chrom, callset, "HWE",    ["EUR", "AFR", "AMR", "SAS", "EAS"], max_idx, min_idx, 1e-6)
    if dataset in ["HGDP", "ABOR", "VANU", "SGDP", "AYTA", "IGDP"]:
        exc_filt = np.full(acc_filt.shape, True)
        hwe_filt = np.full(acc_filt.shape, True)
    #A.6. All filters
    all_filt = acc_filt*pas_filt*snp_filt*max_min_filt*exc_filt*hwe_filt
    #B. CREATE THE ANCESTRAL AND DERIVED ALLELES MAP
    anc_filt, max_min_anc_filt, anc_der_map, ref_alt_alleles = get_anc_der_map(chrom, callset, all_filt, max_min_idx)
    anc_der_map_correct                                      = get_anc_der_map_correct(chrom, callset, ingroup_idx, all_filt, max_min_anc_filt, anc_der_map)
    if vervose:
        print("Filters for chr{} for {} data".format(chrom, dataset), flush = True)
        for x, y in zip([acc_filt, pas_filt, snp_filt, hwe_filt, exc_filt, max_min_filt, anc_filt, max_min_anc_filt, all_filt], ["acc_filt", "pas_filt", "snp_filt", "hwe_filt", "exc_filt", "max_min_filt", "anc_filt", "max_min_anc_filt", "all_filt"]):
            print("{:18}\tshape = {}\tTrue = {}\tperc = {}%".format(y, x.shape, x.sum(), round(100*x.sum()/x.shape[0], 2)), flush=True)

    return all_filt, max_min_anc_filt, anc_der_map, ref_alt_alleles, anc_der_map_correct

def get_outgroup_poly_sites_filt(outgroup_poly_sites_file, chrom, callset, all_filt, max_min_anc_filt):
    out_pol     = allel.SortedIndex(np.unique(np.loadtxt(outgroup_poly_sites_file, usecols=[2], dtype='int')))
    #out_pol    = allel.SortedIndex(          np.loadtxt(outgroup_poly_sites_file, dtype='int'))
    pos         = allel.SortedIndex(callset["{}/variants/POS".format(chrom)])[all_filt][max_min_anc_filt]
    out_filt, _ = pos.locate_intersection(out_pol)
    return ~out_filt