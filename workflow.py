#A. Importing
from templates import *
from variables import *
import itertools
gwf = Workflow()

dataset   = "HGDP"
ancestry  = "Neanderthal"
mean_prob = 0.5
anc1      = "Neanderthal"
anc2      = "Neanderthal"
data1     = "HGDP"
data2     = "HGDP"

#for chromosome in lchromos:
#    gwf.target_from_template("heatmap_{}_{}".format(dataset, chromosome), plot_heatmap(dataset, chromosome))

#for chromosome in lchromos:
#    gwf.target_from_template("{}_chr{}".format(dataset, chromosome), txt_heatmap(dataset, chromosome))


#for population in lpopulat[:2]:
#    gwf.target_from_template("arc_{}_{}_{}_{}".format(dataset, population, ancestry, mean_prob), arc_pop(dataset, population, ancestry, mean_prob))
    
for population in lpopulat:
    gwf.target_from_template("art_{}".format(population), art_gen(dataset, population, ancestry, mean_prob))

for pop1, pop2 in itertools.combinations(lpopulat, 2):
    gwf.target_from_template("int_{}_{}".format(pop1, pop2), intersect(data1, pop1, anc1, data2, pop2, anc2, mean_prob))

for chromosome in lchromos:
    gwf.target_from_template("div_chr{}".format(chromosome), pull(data1, anc1, data2, anc2, mean_prob, chromosome))
