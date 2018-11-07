from gene_graph_lib.compute_complexity import GenomeGraph
import time
from pickle import load

def get_subgraph(input, organism, reference, window=20, start=None, end=None, depth=-1, tails=5, names_list='all'):

    with open(input, 'rb') as f:
        graph = load(f)

    subgraph, aim_chain = graph.generate_subgraph(start, end, reference=reference, window=window, tails=tails, depth=depth)


    
    nodes = set([])
    edges = []

    for org in subgraph:
        for seq in subgraph[org]:
            for i in range(len(seq) - 1):
                edges.append([graph.genes_decode[seq[i]], graph.genes_decode[seq[i+1]], 1])
                nodes.add(graph.genes_decode[seq[i]])
                nodes.add(graph.genes_decode[seq[i+1]])

    ref_chain = [graph.genes_decode[n] for n in subgraph[reference][0]]
    aim_chain = [graph.genes_decode[n] for n in aim_chain]

    return [nodes, edges, aim_chain, ref_chain]



