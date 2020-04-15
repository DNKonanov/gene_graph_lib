import graphviz
from graphviz import Graph
from pygraphviz import AGraph
import math
import time


def draw_graph(input_file, out, ref_chain=[], aim_chain=[], freq_min=1, draw_all=False):
    ps = AGraph(directed=True)
    ps.graph_attr['rankdir'] = 'LR'
    ps.graph_attr['mode'] = 'hier'


    nodes = delete_nodes(ref_chain, input_file[1], freq_min)

    print('Number of nodes: ', end='')
    if draw_all == False:
        print(len(nodes))

    elif draw_all == True:
        print(len(input_file[1]))
    cluster_main = ps.add_subgraph()
    cluster_main.graph_attr['rank'] = '0'


    for i in range(len(ref_chain)):
        shape = 'circle'
        if ref_chain[i] in aim_chain:
            cluster_main.add_node(ref_chain[i], shape=shape, color='red')

        else:

            cluster_main.add_node(ref_chain[i], shape=shape, color='pink')



    for i in input_file[0]:
        if 'PATH' in i:
            shape='box'

        if i in aim_chain or i in ref_chain:
            continue

        if i in nodes:
            ps.add_node(str(i))
            continue
            
        else:
            if draw_all == True:
                ps.add_node(str(i), color='grey')
                continue






    for i in input_file[1]:
        color='black'
        try:
            if ref_chain.index(i[1]) - ref_chain.index(i[0]) == 1:
                color = 'red'
                ps.add_edge(str(i[0]), str(i[1]), color=color, penwidth=str(math.sqrt(i[2])))
                continue

        except ValueError:
            pass

        if i[2] < freq_min:

            if draw_all == True:
                ps.add_edge(str(i[0]), str(i[1]), color='grey', penwidth=str(math.sqrt(i[2])), constraint='false')

            continue

        elif i[0] in nodes and i[1] in nodes:
            ps.add_edge(str(i[0]), str(i[1]), color=color, penwidth=str(math.sqrt(i[2])))

        elif draw_all == True:

            ps.add_edge(str(i[0]), str(i[1]), color='grey', penwidth=str(math.sqrt(i[2])), constraint='false')

    


    ps.draw(out + '.ps', prog='dot')
    ps.write(out + '.dot')







    
def compute_frequence(edges):
    edges.sort()
    for i in range(len(edges) - 1):
        for j in range( i+1,len(edges)):
            if edges[i][0] == edges[j][0] and edges[i][1] == edges[j][1]:
                edges[i][2] += 1
                edges[j][2] = -len(edges)
            else:
                break
    edges_frequency = []
    for i in edges:
        if i[2] > 0:
            edges_frequency.append(i)
    return edges_frequency







def delete_nodes(ref_chain, bonds, freq_min):

    nodes = []

    nodes = set(ref_chain.copy())

    count = 1

    while (count > 0):
        count = 0
        for i in bonds:
            if i[2] < freq_min:
                continue

            if i[0] in nodes and i[1] in nodes:
                continue

            if (i[0] in nodes or i[1] in nodes):
                nodes.add(i[0])
                nodes.add(i[1])
                count += 1

    return list(nodes)


