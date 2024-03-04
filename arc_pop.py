import pandas as pd
import numpy as np
import sys
from variables import *
import os

sys.path.insert(1, '/home/clsj/GenerationInterval/people/moi/scripts')


#########################################################
#                      Load data                        #
#########################################################

dataset = sys.argv[1]
population = sys.argv[2]
ancestry = sys.argv[3]
mean_prob = sys.argv[4]

metadata = pd.read_table("~/GenerationInterval/people/moi/files/metadata.txt")
metadata = metadata[metadata['dat'] == dataset]

#########################################################
#             Make archaic population file              #
#########################################################

data = metadata[metadata['pop'] == population]
samples = data.sample(n=len(data))
pop_ind = (samples['ind'].tolist())
frames = []

for ind in pop_ind:
    human_id = ind
    ind_bed_in = '~/GenerationInterval/people/moi/sandbox/dec/{dataset}/{human_id}.new.txt'.format(dataset = dataset, human_id = human_id)
    ind_file = pd.read_csv(ind_bed_in, sep = '\t')
    ind_file["end"] += 1000 
    frames.append(ind_file)
    pop_df = pd.concat(frames)
    conditions = [
    (pop_df['ancestry'] == 'Vindija') | (pop_df['ancestry'] == 'AmbigNean') | (pop_df['ancestry'] == 'Altai') | (pop_df['ancestry'] == 'Chagyrskaya'),
    (pop_df['ancestry'] == 'Denisova'),
    (pop_df['ancestry'] == 'nonDAVC'),
    (pop_df['ancestry'] == 'Ambiguous'),
    (pop_df['ancestry'] == 'Vindija') | (pop_df['ancestry'] == 'Ambiguous') | (pop_df['ancestry'] == 'Denisova') | (pop_df['ancestry'] == 'nonDAVC') | (pop_df['ancestry'] == 'AmbigNean') | (pop_df['ancestry'] == 'Altai') | (pop_df['ancestry'] == 'Chagyrskaya')
    ]
    values = ['Neanderthal', 'Denisova', 'NonDAVC', 'Ambiguous', 'All']
    pop_df['group'] = np.select(conditions, values)
    anc_df = pop_df.loc[pop_df['group'] == ancestry]
    value=mean_prob
    #fin_df = anc_df.query("mean_prob > @value")
    fin_df = anc_df.query("mean_prob > {mean_prob}".format(mean_prob = mean_prob))

    fin_df.to_csv("{}/{}/archaicpopulations/arc_{}_{}_{}_{}.txt".format(results,ancestry, dataset, population, ancestry, mean_prob), index = False, header = True, sep = '\t')

#with open(os.path.join(completed, 'arc_{}_{}_{}_{}.DONE'.format(dataset, population, ancestry, mean_prob)), 'w') as fp:
    #pass