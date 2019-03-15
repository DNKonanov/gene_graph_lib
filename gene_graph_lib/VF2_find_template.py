from gene_graph_lib.compute_complexity import GenomeGraph
import numpy as np
from networkx.algorithms import isomorphism
import networkx as nx
from networkx.algorithms.isomorphism.isomorphvf2 import DiGraphMatcher
from collections import OrderedDict

def _filter_by_genomes(templates_candidates, g):

    pairs_genomes_sets = {}
    
    print('Create edges genomes set...')
    for genome in g.list_graph:
        for contig in g.list_graph[genome]:
            c = g.list_graph[genome][contig]

            for i in range(len(c) - 1):
                if (c[i], c[i+1]) not in pairs_genomes_sets:
                    pairs_genomes_sets[(c[i], c[i+1])] = set([genome])
                else:
                    pairs_genomes_sets[(c[i], c[i+1])].add(genome)
    
    real_matches = []
    N_genomes = len(g.list_graph)

    for template in templates_candidates:

        new_template = {}
        genomes_check_list = {}
        for edge in template[0]:
            if edge[2] not in genomes_check_list:
                genomes_check_list[edge[2]] = [(edge[0], edge[1], edge[3])]
            else:
                genomes_check_list[edge[2]].append((edge[0], edge[1], edge[3]))

        for genome in genomes_check_list:
            try:
                intersect_all = pairs_genomes_sets[(genomes_check_list[genome][0][0], genomes_check_list[genome][0][1])]
            except KeyError:
                break
            for edge in genomes_check_list[genome]:
                try:
                    intersect_all = intersect_all.intersection(pairs_genomes_sets[(edge[0], edge[1])])
                except:
                    break
                if len(intersect_all) == 0:
                    break
                if len(pairs_genomes_sets[(edge[0], edge[1])]) < edge[2]*N_genomes:
                    
                    intersect_all = []
                    break

                new_template[g.genes_decode[edge[0]]] = (g.genes_decode[edge[1]], pairs_genomes_sets[(edge[0], edge[1])])

            if len(intersect_all) == 0:
                break

        if len(intersect_all) == 0:
            continue

        real_matches.append(new_template)

    return real_matches



def find_template(template, g):
    #crete template adj matrix
    print('Create template adj matrix...')
    clear_edges = template[0]

    Nodes = []
    tCodes = {}
    tDecodes = {}
    code = 0
    for edge in clear_edges:
        if edge[0] not in Nodes:
            Nodes.append(edge[0])
        if edge[0] not in tCodes:
            tCodes[edge[0]] = code
            tDecodes[code] = edge[0]
            code += 1
        
        if edge[1] not in Nodes:
            Nodes.append(edge[1])
        if edge[1] not in tCodes:
            tCodes[edge[1]] = code
            tDecodes[code] = edge[1]
            code += 1




    template_matrix = [[0 for i in range(len(Nodes))] for j in range(len(Nodes))]

    for edge in clear_edges:
        template_matrix[tCodes[edge[0]]][tCodes[edge[1]]] = 1

    #create graph adj matrix
    print('Create genome graph adj matrix...')
    code = 0
    Codes = {}
    Decodes = {}
    for i in g.dict_graph: 
        if i not in Codes:
            
            Codes[i] = code
            Decodes[code] = i
            code += 1
            
        for j in g.dict_graph[i]:
            if j not in Codes:
                
                Codes[j] = code
                Decodes[code] = j
                code += 1 

    GraphMatrix = [[0 for i in range(len(Codes))] for j in range(len(Codes))]
    for gene in g.dict_graph:
        for gene2 in g.dict_graph[gene]:
            
            try:
                GraphMatrix[Codes[gene]][Codes[gene2]] = 1
            except:
                print(gene, gene2)
                continue

    print('Matching...')
    G1 = nx.from_numpy_matrix(np.array(GraphMatrix), create_using=nx.DiGraph())
    G2 = nx.from_numpy_matrix(np.array(template_matrix), create_using=nx.DiGraph())

    GM = DiGraphMatcher(G1,G2)

    hits = [i for i in GM.subgraph_isomorphisms_iter()]
    print('Done!')
    print()
    print('All hits: ', end=' ')

    for j in range(len(hits)):
        hits[j] = {tDecodes[hits[j][i]]: Decodes[i] for i in hits[j]}

    template_candidates = []
    for hit in hits:
        candidate = []

        for edge in template[0]:
            candidate.append((hit[edge[0]], hit[edge[1]], edge[2], edge[3][0]))

        template_candidates.append([candidate, []])

    print(template_candidates)

    filtered_templates = _filter_by_genomes(template_candidates, g)
    return filtered_templates
    

    

    

#A = [[('1','2', 'g0', 0.5), ('2','3', 'g0', 0.3), ('3','4', 'g0', 0.5), ('2', '6', 'g1', 0.1), ('6', '3', 'g1', 0.1)], []]

#g = GenomeGraph()
#g.read_graph('Streptococcus_pneumoniae.sif')
#g.read_graph('../../Escherichia_coli/Escherichia_coli.sif')


#find_template(A, g)