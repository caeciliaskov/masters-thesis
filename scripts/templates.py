#A. Importing
from gwf import Workflow
from variables import *
import itertools

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
    outputs = ['{tables}/div_whole_genome_chr{chromosome}.txt'.format(chromosome = chromosome,tables=tables)]
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
    outputs = ['{}/arc_{}_{}_{}_{}.DONE'.format(completed, dataset, population, ancestry, mean_prob)]
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
    outputs = ['{}/art_{}_{}_{}_{}.DONE'.format(completed, dataset, population, ancestry, mean_prob)]
    options = {'memory'  : '1g', 
               'walltime': '02:00:00', 
               'account' : 'GenerationInterval'}
    spec    = '''
    
    echo “JOBID            : ” $PBS_JOBID
    echo “HOSTNAME         : ” $HOSTNAME
    echo “CONDA_DEFAULT_ENV: ” $CONDA_DEFAULT_ENV 
    
    python arc_pop.py {dataset} {population} {ancestry} {mean_prob}

    rm -f {artificialgenomes}/art_{dataset}_{population}_{ancestry}_{mean_prob}.txt
    touch {artificialgenomes}/art_{dataset}_{population}_{ancestry}_{mean_prob}.txt

    awk '{{if(NR > 1){{print $1"\t"$2"\t"$3"\t"$5"\t"$12}}}}' {archaicpopulations}/arc_{dataset}_{population}_{ancestry}_{mean_prob}.txt \
    | shuf \
    | while read line 
    do 
      bedtools intersect -v -a <(echo "" | awk '{{print "'"${{line}}"'"}}') -b {artificialgenomes}/art_{dataset}_{population}_{ancestry}_{mean_prob}.txt > {artificialgenomes}/art_tmp_{dataset}_{population}_{ancestry}_{mean_prob}.txt
              
      cat {artificialgenomes}/art_tmp_{dataset}_{population}_{ancestry}_{mean_prob}.txt >> {artificialgenomes}/art_{dataset}_{population}_{ancestry}_{mean_prob}.txt
      
      bedtools sort -i {artificialgenomes}/art_{dataset}_{population}_{ancestry}_{mean_prob}.txt > {artificialgenomes}/art_tmp_{dataset}_{population}_{ancestry}_{mean_prob}.txt
      
      mv {artificialgenomes}/art_tmp_{dataset}_{population}_{ancestry}_{mean_prob}.txt {artificialgenomes}/art_{dataset}_{population}_{ancestry}_{mean_prob}.txt
      
    done

    echo ""
    
    touch {completed}/art_{dataset}_{population}_{ancestry}_{mean_prob}.DONE

    '''.format(dataset = dataset, population = population, ancestry = ancestry, mean_prob = mean_prob, artificialgenomes = artificialgenomes, archaicpopulations = archaicpopulations, completed = completed)
    return inputs, outputs, options, spec

def intersect(data1, pop1, anc1, data2, pop2, anc2, mean_prob):
    '''
    A dummy job to just flag that the same job for different inputs have been done
    '''
    inputs  = ['{}/art_{}_{}_{}_{}.DONE'.format(completed, data1, pop1, anc1, mean_prob), '{}/art_{}_{}_{}_{}.DONE'.format(completed, data2, pop2, anc2, mean_prob)]
    outputs = ['{}/int_{}_{}_{}_{}_{}_{}_{}.DONE'.format(completed, data1, pop1, anc1, data2, pop2, anc2, mean_prob)]
    options = {'memory'  : '1g', 
               'walltime': '01:00:00', 
               'account' : 'GenerationInterval'}
    spec    = '''
    
    echo “JOBID            : ” $PBS_JOBID
    echo “HOSTNAME         : ” $HOSTNAME
    echo “CONDA_DEFAULT_ENV: ” $CONDA_DEFAULT_ENV 
    
    bedtools intersect -wa -wb -a {artificialgenomes}/art_{data1}_{pop1}_{anc1}_{mean_prob}.txt -b {artificialgenomes}/art_{data2}_{pop2}_{anc2}_{mean_prob}.txt > {intersection}/ind_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt
    
    bedtools intersect -a <(awk '{{print $1"\t"$2"\t"$3}}' {artificialgenomes}/art_{data1}_{pop1}_{anc1}_{mean_prob}.txt) -b <(awk '{{print $1"\t"$2"\t"$3}}' {artificialgenomes}/art_{data2}_{pop2}_{anc2}_{mean_prob}.txt) > {intersection}/int_tmp_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt
    
    paste {intersection}/int_tmp_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt <(awk '{{print $5"\t"$10}}' {intersection}/ind_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt) > {intersection}/int_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt
    
    bedtools sort -i {intersection}/int_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt > {intersection}/int_tmp_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt
    
    cat {intersection}/int_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt | tr ' ' '\t' > {intersection}/int_tmp_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt 
    
    cat {intersection}/int_tmp_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt > {intersection}/int_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt

    rm {intersection}/ind_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt {intersection}/int_tmp_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.txt

    echo ""

    touch {completed}/int_{data1}_{pop1}_{anc1}_{data2}_{pop2}_{anc2}_{mean_prob}.DONE

    '''.format(data1 = data1, pop1 = pop1, anc1 = anc1, data2 = data2, pop2 = pop2, anc2 = anc2, mean_prob = mean_prob, artificialgenomes = artificialgenomes, intersection = intersection, completed = completed)
    return inputs, outputs, options, spec

def pull(data1, anc1, data2, anc2, mean_prob, chrom):
    '''
    A dummy job to just flag that the same job for different inputs have been done
    '''
    inputs = []
    for pop1, pop2 in itertools.combinations(lpopulat, 2):
        inputs.append('{}/int_{}_{}_{}_{}_{}_{}_{}.DONE'.format(completed, data1,pop1,anc1,data2,pop2,anc2,mean_prob))
    outputs = ['{}/div_{}_{}_{}_{}_{}_chr{}.DONE'.format(completed, data1, anc1, data2, anc2, mean_prob, chrom)]
    options = {'memory'  : '50g', 
               'walltime': '20:00:00', 
               'account' : 'GenerationInterval'}
    spec    = '''

    echo “JOBID            : ” $PBS_JOBID
    echo “HOSTNAME         : ” $HOSTNAME
    echo “CONDA_DEFAULT_ENV: ” $CONDA_DEFAULT_ENV 

    python pullvariants.py {} {} {} {} {} {}

    '''.format(data1, anc1, data2, anc2, mean_prob, chrom)
    return inputs, outputs, options, spec
