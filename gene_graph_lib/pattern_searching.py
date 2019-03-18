import networkx
import collections
from gene_graph_lib import compute_complexity
import time
from gene_graph_lib.compute_complexity import GenomeGraph
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
