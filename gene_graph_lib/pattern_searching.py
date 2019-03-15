import networkx
import collections
from gene_graph_lib import compute_complexity
import time
from gene_graph_lib.compute_complexity import GenomeGraph
# import matplotlib.pyplot as plt
# from compute_complexity import GenomeGraph




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


class SpyderPattern(object):
    def __init__(
            self,
            min_legs_len=50,
            min_legs_num=4
    ):
        self.name = 'Spider'
        self.min_legs_len = min_legs_len,
        self.min_legs_num = min_legs_num


class PenguinPatter(object):
    def __init__(
            self,
            min_distance=1,
            min_length=1,
            max_length=20
    ):
        self.name = 'Penguin',
        self.min_distance = min_distance,
        self.max_length = max_length,
        self.min_length = min_length,


def count_garlic_for_gen(
        genome_dict,
        genome,
        hashed,
        garlic_pattern) -> list:

    def garlic_number(
            left_contig,
            right_contig,
            min_insert_per_pos,
            max_insert_per_pos) -> list:
        _result = []
        for left_gen, right_gen in zip(left_contig[:-1], left_contig[1:]):
            if (left_gen, right_gen) in hashed:
                _result.append((left_gen, right_gen))
                continue
            if (left_gen not in right_contig) or (
                    right_gen not in right_contig):
                continue
            dif = right_contig.index(right_gen) - right_contig.index(left_gen)
            if (dif > min_insert_per_pos) and (dif < max_insert_per_pos):
                _result.append((left_gen, right_gen))
                hashed.add((left_gen, right_gen))
        return _result

    result = []
    for contig_name in genome_dict[genome]:
        current_contig = genome_dict[genome][contig_name]
        for other_genome in genome_dict.keys() - set(genome):
            for other_contig in genome_dict[other_genome].values():
                result.extend(garlic_number(
                    current_contig,
                    other_contig,
                    min_insert_per_pos=garlic_pattern.min_insert_per_pos,
                    max_insert_per_pos=garlic_pattern.max_insert_per_pos))
    return result


def count_spyder_for_gen(
        genome_dict,
        genome,
        hashed,
        spider_pattern) -> list:
    pass


def count_pengiun_for_gen(
        genome_dict,
        genome,
        hashed,
        penguin_pattern) -> list:
    pass


pattern_searchers = {
    'Garlic': count_garlic_for_gen,
}


def count_patten_in_gen_for_each(
        genome_dict,
        pattern) -> dict:
    result = {}
    hashed = set()

    count_pattern_for_gen = pattern_searchers[pattern.name]

    for genome, contigs in genome_dict.items():
        result[genome] = count_pattern_for_gen(
            genome_dict,
            genome,
            hashed,
            pattern
        )

    for k, v in result.items():
        result[k] = {
            pr: wght for pr, wght in dict(collections.Counter(v)).items()
            if (wght >= pattern.min_insert_count) and
               (wght <= pattern.max_insert_per_pos)
        }
    return result


def find_garlic(g, min_distinct_paths=5, min_in_weight=0.2, min_out_weight=0.4, max_insert_len=25, min_insert_len=5):

    hits = {}
    in_weight = int(min_in_weight*len(g.list_graph))
    out_weight = int(min_out_weight*len(g.list_graph))

    for genome in g.list_graph:
        print(genome)
        for contig in g.list_graph[genome]:
            c = g.list_graph[genome][contig]

            for i in range(1, len(c) - 2):
                start, end = c[i], c[i+1]
                
                try:
                    if g.edges_weights[(c[i-1], c[i])] < out_weight or g.edges_weights[(c[i+1], c[i+2])] < out_weight:
                        continue

                    if g.edges_weights[(c[i], c[i+1])] < in_weight:
                        hits[(start, end)] = set([])
                        continue
                
                except:
                    continue

                if (start, end) in hits:
                    continue

                print(i, end='')
                if (start, end) not in hits:
                    hits[(start, end)] = set([])

                for other_genome in g.list_graph:
                    for other_contig in g.list_graph[other_genome]:
                        other_c = g.list_graph[other_genome][other_contig]

                        if start not in other_c or end not in other_c:
                            continue

                        left_index, rigth_index = other_c.index(start), other_c.index(end)
                        if left_index > rigth_index:
                            left_index, rigth_index = rigth_index, left_index
                        
                        path = tuple(other_c[left_index:rigth_index + 1])

                        if len(path) > max_insert_len or len(path) < min_insert_len:
                            continue

                        if path in hits[(start, end)]:
                            continue
                        
                        hits[(start, end)].add(path)
                print('\r', end='')


    hits = [(g.genes_decode[h[0]], g.genes_decode[h[1]]) for h in hits if len(hits[h]) >= min_distinct_paths]

    return hits



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
    
    all_penguins = {}
    for h in hits:
        for e in h[:-1]:
            if e not in all_penguins:
                all_penguins[e] = 1
            else:
                all_penguins[e] += 1
    return all_penguins

def find_spyder(g, min_legs_len=50, min_legs_num=4):
    

    all_spyders = {}
    no_spyders = set([])

    for genome in g.list_graph:
        print(genome)
        for contig in g.list_graph[genome]:
            c = g.list_graph[genome][contig]

            for i in range(len(c)):

                if c[i] in all_spyders or c[i] in no_spyders:
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
                    all_spyders[c[i]] = len(uniq)

                else:
                    no_spyders.add(c[i])

    all_spyders = {g.genes_decode[i]: all_spyders[i] for i in all_spyders}
    return all_spyders