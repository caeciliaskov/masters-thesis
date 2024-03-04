#A. Importing
from gwf import Workflow
from variables import *
import numpy as np
import vcf_processing_functions as vpf
import allel
import zarr
import pandas as pd

#B. Templates
##B.1.
def plot_heatmap(dataset, chromosome):
    '''
    A dummy job to just flag that the same job for different inputs have been done
    '''
    inputs  = []
    outputs = ['chr{chromosome}.png'.format(chromosome = chromosome)]
    options = {'memory'  : '1g', 
               'walltime': '01:00:00', 
               'account' : 'GenerationInterval'}
    spec    = '''

    echo “JOBID            : ” $PBS_JOBID
    echo “HOSTNAME         : ” $HOSTNAME
    echo “CONDA_DEFAULT_ENV: ” $CONDA_DEFAULT_ENV 

    python divergence.py {dataset} {chromosome}

    '''.format(dataset = dataset, chromosome = chromosome)
    return inputs, outputs, options, spec


def txt_heatmap(dataset, chromosome):
    '''
    A dummy job to just flag that the same job for different inputs have been done
    '''
    inputs  = []
    outputs = ['{results}/tables/div_Human_chr{chromosome}.txt'.format(chromosome = chromosome,results=results)]
    options = {'memory'  : '5g', 
               'walltime': '05:00:00', 
               'account' : 'GenerationInterval'}
    spec    = '''

    echo “JOBID            : ” $PBS_JOBID
    echo “HOSTNAME         : ” $HOSTNAME
    echo “CONDA_DEFAULT_ENV: ” $CONDA_DEFAULT_ENV 

    python divergence_txt.py {dataset} {chromosome}

    '''.format(dataset = dataset, chromosome = chromosome)
    return inputs, outputs, options, spec

def arc_pop(dataset, population, ancestry, mean_prob):
    '''
    A dummy job to just flag that the same job for different inputs have been done
    '''
    inputs  = []
    outputs = ['{}/{}/completed/arc_{}_{}_{}_{}.DONE'.format(results,ancestry, dataset, population, ancestry, mean_prob)]
    options = {'memory'  : '1g', 
               'walltime': '00:10:00', 
               'account' : 'GenerationInterval'}
    spec    = '''

    echo “JOBID            : ” $PBS_JOBID
    echo “HOSTNAME         : ” $HOSTNAME
    echo “CONDA_DEFAULT_ENV: ” $CONDA_DEFAULT_ENV 

    python arc_pop.py {dataset} {population} {ancestry} {mean_prob}

    '''.format(dataset = dataset, population = population, ancestry = ancestry, mean_prob = mean_prob)
    return inputs, outputs, options, spec

def art_gen(dataset, population, ancestry, mean_prob):
    '''
    A dummy job to just flag that the same job for different inputs have been done
    '''
    inputs  = []
    outputs = ['{}/{}/completed/art_{}_{}_{}_{}.DONE'.format(results,ancestry, dataset, population, ancestry, mean_prob)]
    options = {'memory'  : '1g', 
               'walltime': '02:00:00', 
               'account' : 'GenerationInterval'}
    spec    = '''
    
    echo “JOBID            : ” $PBS_JOBID
    echo “HOSTNAME         : ” $HOSTNAME
    echo “CONDA_DEFAULT_ENV: ” $CONDA_DEFAULT_ENV 
    
    python arc_pop.py {dataset} {population} {ancestry} {mean_prob}

    rm -f {results}/{ancestry}/artificialgenomes/art_{dataset}_{population}_{ancestry}_{mean_prob}.txt
    touch {results}/{ancestry}/artificialgenomes/art_{dataset}_{population}_{ancestry}_{mean_prob}.txt

    awk '{{if(NR > 1){{print $1"\t"$2"\t"$3"\t"$5"\t"$12}}}}' {results}/{ancestry}/archaicpopulations/arc_{dataset}_{population}_{ancestry}_{mean_prob}.txt \
    | shuf \
    | while read line 
    do 
      bedtools intersect -v -a <(echo "" | awk '{{print "'"${{line}}"'"}}') -b {results}/{ancestry}/artificialgenomes/art_{dataset}_{population}_{ancestry}_{mean_prob}.txt > {results}/{ancestry}/artificialgenomes/art_tmp_{dataset}_{population}_{ancestry}_{mean_prob}.txt
              
      cat {results}/{ancestry}/artificialgenomes/art_tmp_{dataset}_{population}_{ancestry}_{mean_prob}.txt >> {results}/{ancestry}/artificialgenomes/art_{dataset}_{population}_{ancestry}_{mean_prob}.txt
      
      bedtools sort -i {results}/{ancestry}/artificialgenomes/art_{dataset}_{population}_{ancestry}_{mean_prob}.txt > {results}/{ancestry}/artificialgenomes/art_tmp_{dataset}_{population}_{ancestry}_{mean_prob}.txt
      
      mv {results}/{ancestry}/artificialgenomes/art_tmp_{dataset}_{population}_{ancestry}_{mean_prob}.txt {results}/{ancestry}/artificialgenomes/art_{dataset}_{population}_{ancestry}_{mean_prob}.txt
      
    done

    echo ""
    
    touch {results}/{ancestry}/completed/art_{dataset}_{population}_{ancestry}_{mean_prob}.DONE

    '''.format(dataset = dataset, population = population, ancestry = ancestry, mean_prob = mean_prob, results = results)
    return inputs, outputs, options, spec

def intersect(data1, pop1, anc1, data2, pop2, anc2, mean_prob):
    '''
    A dummy job to just flag that the same job for different inputs have been done
    '''
    inputs  = ['{}/{}/completed/art_{}_{}_{}_{}.DONE'.format(results,anc1, data1, pop1, anc1, mean_prob), '{}/{}/completed/art_{}_{}_{}_{}.DONE'.format(results,anc1, data2, pop2, anc2, mean_prob)]
    outputs = ['{}/{}/completed/int_{}_{}_{}_{}_{}_{}_{}.DONE'.format(results,anc1, data1, pop1, anc1, data2, pop2, anc2, mean_prob)]
    options = {'memory'  : '1g', 
               'walltime': '01:00:00', 
               'account' : 'GenerationInterval'}
    spec    = '''
    
    echo “JOBID            : ” $PBS_JOBID
    echo “HOSTNAME         : ” $HOSTNAME
    echo “CONDA_DEFAULT_ENV: ” $CONDA_DEFAULT_ENV 
    
    bedtools intersect -wa -wb -a {results}/{ancestry}/artificialgenomes/art_{data1}_{pop1}_{anc1}_{mean_prob}.txt -b {results}/{ancestry}/artificialgenomes/art_{data2}_{pop2}_{anc2}_{mean_prob}.txt > {results}/{ancestry}/intersection/ind_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt
    
    bedtools intersect -a <(awk '{{print $1"\t"$2"\t"$3}}' {results}/{ancestry}/artificialgenomes/art_{data1}_{pop1}_{anc1}_{mean_prob}.txt) -b <(awk '{{print $1"\t"$2"\t"$3}}' {results}/{ancestry}/artificialgenomes/art_{data2}_{pop2}_{anc2}_{mean_prob}.txt) > {results}/{ancestry}/intersection/int_tmp_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt
    
    paste {results}/{ancestry}/intersection/int_tmp_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt <(awk '{{print $5"\t"$10}}' {results}/{ancestry}/intersection/ind_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt) > {results}/{ancestry}/intersection/int_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt
    
    bedtools sort -i {results}/{ancestry}/intersection/int_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt > {results}/{ancestry}/intersection/int_tmp_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt
    
    cat {results}/{ancestry}/intersection/int_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt | tr ' ' '\t' > {results}/{ancestry}/intersection/int_tmp_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt 
    
    cat {results}/{ancestry}/intersection/int_tmp_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt > {results}/{ancestry}/intersection/int_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt

    rm {results}/{ancestry}/intersection/ind_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt {results}/{ancestry}/intersection/int_tmp_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt

    echo ""

    touch {results}/{ancestry}/completed/int_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.DONE

    '''.format(data1 = data1, pop1 = pop1, anc1 = anc1, data2 = data2, pop2 = pop2, anc2 = anc2, mean_prob = mean_prob, results=results,ancestry=anc1)
    return inputs, outputs, options, spec

def pull(data1, anc1, data2, anc2, mean_prob, chrom):
    '''
    A dummy job to just flag that the same job for different inputs have been done
    '''
    inputs = []
    for i in range(len(lpopulat)):
        for j in range(i, len(lpopulat)):
            pop1 = lpopulat[i]
            pop2 = lpopulat[j]
            inputs.append('{}/{}/completed/int_{}_{}_{}_{}_{}_{}_{}.DONE'.format(results,anc1, data1, pop1, anc1, data2, pop2, anc2, mean_prob))
    outputs = []
    for i in range(len(lpopulat)):
        for j in range(i, len(lpopulat)):
            pop1 = lpopulat[i]
            pop2 = lpopulat[j]
            outputs.append('{}/{}/completed/div_{}_{}_{}_{}_{}_{}_{}_chr{}.DONE'.format(results,anc1, data1, pop1, anc1, data2, pop2, anc2, mean_prob,chrom))
    options = {'memory'  : '30g', 
               'walltime': '10:00:00', 
               'account' : 'GenerationInterval'}
    spec    = '''

    echo “JOBID            : ” $PBS_JOBID
    echo “HOSTNAME         : ” $HOSTNAME
    echo “CONDA_DEFAULT_ENV: ” $CONDA_DEFAULT_ENV 

    python pullvariants_MOI2.py {} {} {} {} {} {}

    '''.format(data1, anc1, data2, anc2, mean_prob, chrom)
    return inputs, outputs, options, spec

def genome_length():
    '''
    A dummy job to just flag that the same job for different inputs have been done
    '''
    inputs = []
    outputs = ['genome_length_Human.txt']
    options = {'memory'  : '5g', 
               'walltime': '01:00:00', 
               'account' : 'GenerationInterval'}
    spec    = '''

    echo “JOBID            : ” $PBS_JOBID
    echo “HOSTNAME         : ” $HOSTNAME
    echo “CONDA_DEFAULT_ENV: ” $CONDA_DEFAULT_ENV 

    python genome_length.py

    '''
    return inputs, outputs, options, spec

def all_done(inputs, output):
        '''
        A dummy job to just flag that the same job for different inputs have been done
        '''
        outputs = [output]
        options = {'memory' : '1g', 'walltime': '00:01:00', 'account' : 'GenerationInterval'}
        spec    = '''touch {output}'''.format(output = output)
        return inputs, outputs, options, spec