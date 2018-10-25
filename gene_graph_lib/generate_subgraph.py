from gene_graph_lib.compute_complexity import GenomeGraph
import time
from pickle import load

def get_subgraph(input, organism, reference_stamm, window=20, start=None, end=None, depth=-1, tails=5, names_list='all'):

    with open(input, 'rb') as f:
        graph = load(f)

    subgraph, aim_chain = graph.generate_subgraph(start, end, reference=reference_stamm, window=window, tails=tails, depth=depth)

    f_out = ''
    f_freq = ''

    
    in_aim = 0
    is_ref = 0
    for stamm in subgraph:
        for contig in subgraph[stamm]:
            for i in range(len(contig) - 1):

                try:
                    if aim_chain.index(contig[i + 1]) - aim_chain.index(contig[i]) == 1:
                        in_aim = 1
                except ValueError:
                    pass

                try:
                    if subgraph[reference_stamm][0].index(contig[i + 1]) - subgraph[reference_stamm][0].index(contig[i]) == 1:
                        is_ref = 1
                except ValueError:
                    pass


                line = graph.genes_decode[contig[i]] + ' ' + graph.genes_decode[contig[i + 1]] + ' ' + stamm + ' ' + str(in_aim) + ' ' + str(is_ref) + '\n'
                f_out += line

                in_aim = 0
                is_ref = 0

    freq = {}

    for name in subgraph:
        for contig in subgraph[name]:
            for i in range(len(contig) - 1):
                if graph.genes_decode[contig[i]] + ' ' + graph.genes_decode[contig[i + 1]] not in freq:
                    freq[graph.genes_decode[contig[i]] + ' ' + graph.genes_decode[contig[i + 1]]] = [name]

                else:
                    freq[graph.genes_decode[contig[i]] + ' ' + graph.genes_decode[contig[i + 1]]].append(name)


    for pair in freq:
        f_freq += pair + ' ' + str(len(freq[pair])) + ' '
        for name in freq[pair]:
            f_freq += name
            if name == freq[pair][-1]:
                f_freq += '\n'
                continue
            f_freq +='|'

    return f_out, f_freq



