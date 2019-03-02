import networkx
import collections
from gene_graph_lib import compute_complexity
import time

import matplotlib.pyplot as plt
from compute_complexity import GenomeGraph

class GarlicPattern(object):
    def __init__(
            self,
            min_insert_per_pos=1,
            max_insert_per_pos=11,
            min_insert_count=1,
            max_insert_count=10,
            left_weight=.25,
            midl_weight=.25,
            right_weight=.25
    ):
        self.name = 'Garlic'
        self.min_insert_per_pos = min_insert_per_pos
        self.max_insert_per_pos = max_insert_per_pos
        self.min_insert_count = min_insert_count
        self.max_insert_count = max_insert_count
        self.left_weight = left_weight
        self.midl_weight = midl_weight
        self.right_weight = right_weight


def garlic_number(
        left_contig,
        right_contig,
        min_insert_per_pos,
        max_insert_per_pos,
        hashed) -> list:
    result = []
    for left_gen, right_gen in zip(left_contig[:-1], left_contig[1:]):
        if (left_gen, right_gen) in hashed:
            result.append((left_gen, right_gen))
            continue
        if (left_gen not in right_contig) or (right_gen not in right_contig):
            continue
        dif = right_contig.index(right_gen) - right_contig.index(left_gen)
        if (dif > min_insert_per_pos) and (dif < max_insert_per_pos):
            result.append((left_gen, right_gen))
            hashed.add((left_gen, right_gen))
    return result


def count_garlic_for_gen(
        genome_dict,
        genome,
        hashed,
        garlic_pattern) -> list:
    result = []
    for contig_name in genome_dict[genome]:
        current_contig = genome_dict[genome][contig_name]
        for other_genome in genome_dict.keys() - set(genome):
            for other_contig in genome_dict[other_genome].values():
                result.extend(garlic_number(
                    current_contig,
                    other_contig,
                    min_insert_per_pos=garlic_pattern.min_insert_per_pos,
                    max_insert_per_pos=garlic_pattern.max_insert_per_pos,
                    hashed=hashed))
    return result


def count_garlic_in_gen_for_each(
        genome_dict,
        garlic_pattern=GarlicPattern()) -> dict:
    result = {}
    hashed = set()
    for genome, contigs in genome_dict.items():
        result[genome] = count_garlic_for_gen(
            genome_dict,
            genome,
            hashed,
            garlic_pattern=garlic_pattern
        )
    for k, v in result.items():
        result[k] = {
            pr: wght for pr, wght in dict(collections.Counter(v)).items()
            if (wght >= garlic_pattern.min_insert_count) and
               (wght <= garlic_pattern.max_insert_per_pos)
        }

    return result


# test_genome_sample = {
#     'a': {1: [1, 2, 3, 4, 5]},
#     'b': {1: [1, 3, 5]},
#     'c': {1: [1, 3, 4, 5]}
# }
#
# print(count_garlic_in_gene(test_genome_sample))

#test_graph = compute_complexity.GenomeGraph()
#test_graph.read_graph('../Streptococcus_pneumoniae.sif')
#print('readed')
# for k, v in test_graph.list_graph.items():
#     print(k, list(map(len, v.values())))

#keys_ = [k for k in test_graph.list_graph.keys()][:15]

#t = time.time()
#for k, v in count_garlic_in_gene({k: test_graph.list_graph[k] for k in keys_}).items():
#    print(k, len(v))
#print(time.time() - t)






def find_transposition(g, max_length=20, min_length=1, min_distance=5, max_distance=100000, conservativity=0.5):
    pairs_genomes_sets = {}

    cons_limit = conservativity*len(g.list_graph)


    hits = {}
    
    print('Create edges genomes set...')
    for genome in g.list_graph:
        for contig in g.list_graph[genome]:
            c = g.list_graph[genome][contig]

            for i in range(len(c) - 1):
                if (c[i], c[i+1]) not in pairs_genomes_sets:
                    pairs_genomes_sets[(c[i], c[i+1])] = set([genome])
                else:
                    pairs_genomes_sets[(c[i], c[i+1])].add(genome)


    for genome in g.list_graph:
        print(genome)
        for contig in g.list_graph[genome]:
            c = g.list_graph[genome][contig]

            for i in range(len(c) - 1):
                if len(g.dict_graph[c[i]]) > 1:
                    for gene in g.dict_graph[c[i]]:
                        if gene == c[i+1]:
                            continue

                        try:
                            if (min_distance <= abs(c.index(gene) - i) <= max_distance) == False:
                                continue

                            tr_index = c.index(gene)
                            
                            for j in range(tr_index-max_length, tr_index + max_length + 1):
                                if (c[j] in g.dict_graph[c[tr_index]]) == False:
                                    continue

                                if (c[i+1] in g.dict_graph[c[j]]) == False:
                                    continue

                                if abs(i - j) < min_distance:
                                    continue
                                
                                if(abs(j-tr_index)) < min_length:
                                    continue

                                intersect_all = pairs_genomes_sets[(c[i], c[i+1])]

                                for k in range(min(j, tr_index), max(tr_index, j)):
                                    intersect_all = intersect_all.intersection(pairs_genomes_sets[(c[k], c[k+1])])
                                
                                
                                if len(intersect_all) >= cons_limit:
                                    
                                    if genome not in hits:
                                        hits[genome] = set([(
                                            g.genes_decode[c[i]], 
                                            g.genes_decode[c[i+1]], 
                                            g.genes_decode[c[j]], 
                                            g.genes_decode[c[tr_index]],
                                            abs(j-tr_index))])

                                    else:
                                        hits[genome].add((
                                            g.genes_decode[c[i]], 
                                            g.genes_decode[c[i+1]], 
                                            g.genes_decode[c[j]], 
                                            g.genes_decode[c[tr_index]],
                                            abs(j-tr_index)))

                        except:
                            continue
    return hits



#test_graph = GenomeGraph()
#test_graph.read_graph('../../data/Streptococcus_pneumoniae/Streptococcus_pneumoniae.sif')
#test_graph.read_graph('../../data/Escherichia_coli/Escherichia_coli.sif')
#test_graph.read_graph('../../data/Streptococcus_pneumoniae/Streptococcus_pneumoniae.sif')
#test_graph.read_graph('../../data/Neisseria_gonorrhoeae.sif')


#E = find_transposition(test_graph, conservativity=0.0)


#count = 0

#for genome in E:

 #   for contig in test_graph.list_graph[genome]:

#        X = {gene: 0 for gene in test_graph.list_graph[genome][contig]}
#        print(X)
#        for e in E[genome]:
#            try:
#                X[test_graph.genes_code[e[0]]] += 1
#                X[test_graph.genes_code[e[1]]] += 1
#            except:
#                continue

    

#    if count > 20:
#        break


#plt.plot([X[i] for i in X])
#plt.show()

#count += 1


#U = []
#for genome in E:
#    U += [i[4] for i in E[genome]]

#plt.hist(U)
#plt.show()





def find_spider(g, min_legs_len=50, min_legs_num=4):
    

    all_spiders = set([])
    no_spiders = set([])

    for genome in g.list_graph:
        print(genome)
        for contig in g.list_graph[genome]:
            c = g.list_graph[genome][contig]

            for i in range(len(c)):

                if c[i] in all_spiders or c[i] in no_spiders:
                    continue

                uniq = set([])
                if len(g.dict_graph) < min_legs_num:
                    continue
                for other_genome in g.list_graph:
                    for other_contig in g.list_graph[other_genome]:
                        other_c = g.list_graph[other_genome][other_contig]

                        if c[i] not in other_c:
                            continue
                        k = other_c.index(c[i])

                        if len(set(other_c[k-min_legs_len:k+min_legs_len+1]).intersection(set(c[i-min_legs_len:i+min_legs_len+1]))) == 1:
                            uniq.add(tuple(other_c[k-min_legs_len:k]))

                if len(uniq) >= min_legs_num:
                    all_spiders.add(c[i])

                else:
                    no_spiders.add(c[i])

    print(len(all_spiders))
    all_spiders = set([g.genes_decode[i] for i in all_spiders])
    for a in all_spiders:
        print(a)
