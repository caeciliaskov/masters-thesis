#A. Importing
from templates import *
from variables import *
gwf = Workflow()

dataset   = "HGDP"
mean_prob = 0.5
data1     = "HGDP"
data2     = "HGDP"

#for chromosome in lchromos:
#    gwf.target_from_template("heatmap_{}_{}".format(dataset, chromosome), plot_heatmap(dataset, chromosome))

for chromosome in lchromos:
    gwf.target_from_template("d_Human_{}".format(chromosome), txt_heatmap(dataset, chromosome))

#for population in lpopulat:
#    gwf.target_from_template("arc_{}_{}_{}_{}".format(dataset, population, ancestry, mean_prob), arc_pop(dataset, population, ancestry, mean_prob))

for ancestry in ["Neanderthal","Denisova"]:
    anc1 = ancestry
    anc2 = ancestry
    inputs = []
    for population in lpopulat:
        gwf.target_from_template("a_{}_{}".format(ancestry, population), art_gen(dataset, population, ancestry, mean_prob))
        inputs.append('{}/{}/completed/art_{}_{}_{}_{}.DONE'.format(results,ancestry, dataset, population, ancestry, mean_prob))
    gwf.target_from_template("a_{}".format(ancestry), all_done(inputs = inputs, output = '{}/completed/art_{}.DONE'.format(results, ancestry)))

    inputs = []
    for i in range(len(lpopulat)):
        for j in range(i, len(lpopulat)):
            pop1 = lpopulat[i]
            pop2 = lpopulat[j]
            gwf.target_from_template("i_{}_{}_{}".format(ancestry, pop1, pop2), intersect(data1, pop1, anc1, data2, pop2, anc2, mean_prob))
            inputs.append('{}/{}/completed/int_{}_{}_{}_{}_{}_{}_{}.DONE'.format(results,ancestry, data1, pop1, anc1, data2, pop2, anc2, mean_prob))
    gwf.target_from_template("i_{}".format(ancestry), all_done(inputs = inputs, output = '{}/completed/int_{}.DONE'.format(results, ancestry)))

    inputs = []
    for chromosome in lchromos:
        gwf.target_from_template("d_{}_{}".format(ancestry, chromosome), pull(data1, anc1, data2, anc2, mean_prob, chromosome))
        for i in range(len(lpopulat)):
            for j in range(i, len(lpopulat)):
                pop1 = lpopulat[i]
                pop2 = lpopulat[j]
                inputs.append('{}/{}/completed/div_{}_{}_{}_{}_{}_{}_{}_chr{}.DONE'.format(results,ancestry, data1, pop1, anc1, data2, pop2, anc2, mean_prob, chromosome))
    gwf.target_from_template("d_{}".format(ancestry), all_done(inputs = inputs, output = '{}/completed/div_{}.DONE'.format(results, ancestry)))

