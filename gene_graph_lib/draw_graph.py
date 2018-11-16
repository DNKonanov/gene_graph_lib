import math
import time
from pygraphviz import AGraph

def get_json_graph(input_subgraph, freq_min, da=False):

    genes = list(input_subgraph[0])
    bonds = compute_frequence(input_subgraph[1])
    aim_chain = input_subgraph[2]
    ref_chain = input_subgraph[3]

    print('Rendering...')
    json_graph = to_json_representation([genes, bonds], ref_chain=ref_chain,
                                        aim_chain=aim_chain, freq_min=freq_min, draw_all=da)

    
    positions, edge_points = get_graphviz_layout([genes, bonds], ref_chain=ref_chain,
                                        aim_chain=aim_chain, freq_min=freq_min, draw_all=da)

                     
    for node in json_graph['nodes']:
        try:
            node['data']['position'] = {'x': float(positions[node['data']['id']].split(',')[0]), 'y': float(positions[node['data']['id']].split(',')[1])}
        except:
            node['data']['position'] = {'x': 0, 'y': 0}
            continue


    return json_graph

# Returns JSON representation of subgraph ready for frontend consumption



def get_graphviz_layout(input_subgraph, ref_chain=[], aim_chain=[], freq_min=1, draw_all=False):
    print('SD')
    ps = AGraph(directed=True)
    ps.graph_attr['rankdir'] = 'LR'
    ps.graph_attr['mode'] = 'hier'


    nodes = delete_nodes(ref_chain, input_subgraph[1], freq_min)

    print('Number of nodes: ', end='')
    if draw_all == False:
        print(len(nodes))

    elif draw_all == True:
        print(len(input_subgraph[1]))
    cluster_main = ps.add_subgraph()
    cluster_main.graph_attr['rank'] = '0'


    for i in range(len(ref_chain)):
        shape = 'circle'
        if ref_chain[i] in aim_chain:
            cluster_main.add_node(ref_chain[i], shape=shape, color='red')

        else:

            cluster_main.add_node(ref_chain[i], shape=shape, color='pink')



    for i in input_subgraph[0]:
        if i in aim_chain or i in ref_chain:
            continue

        if i in nodes:
            ps.add_node(str(i))
            continue
            
        else:
            if draw_all == True:
                ps.add_node(str(i), color='grey')
                continue

    for i in input_subgraph[1]:
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

    ps.layout(prog='dot')

    positions = {n: n.attr['pos'] for n in ps.nodes()}
    edge_points = {edge: edge.attr['pos'] for edge in ps.edges()}
    return positions, edge_points

    

def to_json_representation(input_subgraph, ref_chain=[], aim_chain=[], freq_min=1, draw_all=False):

    print('ND')

    # Deleting nodes that occures less than freq_min times
    nodes = delete_nodes(ref_chain, input_subgraph[1], freq_min)
    max_width = max([edge[2] for edge in input_subgraph[1]])

    print('Number of nodes: ', end='')

    if draw_all == False:
        print(len(nodes))
    elif draw_all == True:
        print(len(input_subgraph[1]))

    nodes_list = []
    if len(aim_chain) > 0:
        nodes_list.append({'data':
                                {'id': 'm', 'color': '#d3d3d3'}})

    for i in range(len(ref_chain)):
        shape = 'circle'
        if ref_chain[i] in aim_chain:
            nodes_list.append({'data':
                {'id': ref_chain[i], 'shape': shape, 'color': '#ff0000', 'parent':'m'}})
        else:
            nodes_list.append({'data':
                {'id': ref_chain[i], 'shape': shape, 'color': 'pink'}})

    for i in input_subgraph[0]:

        if i in aim_chain or i in ref_chain:
            continue

        if i in nodes:
            nodes_list.append({'data':{'id': str(i), 'shape': shape, 'color': 'green'}})
            continue

        else:
            if draw_all == True:
                nodes_list.append({'data':
                    {'id': str(i), 'shape': shape, 'color': 'grey'}})
                continue
    # print(nodes_list)

    edges_list = []
    for i in input_subgraph[1]:
        color = 'black'
        try:
            if ref_chain.index(i[1]) - ref_chain.index(i[0]) == 1:
                color = 'red'
                edges_list.append({'data':{'source': str(i[0]), 'target': str(i[1]),
                                    'color': '#ff0000', 'penwidth': str(10*math.sqrt(i[2]/max_width))}})
                continue
        except ValueError:
            pass

        if i[2] < freq_min:
            if draw_all == True:
                edges_list.append({'data':{'source': str(i[0]), 'target': str(i[1]),
                                   'color': 'grey', 'penwidth': str(10*math.sqrt(i[2]/max_width))}})

            continue

        elif i[0] in nodes and i[1] in nodes:
            edges_list.append({'data':{'source': str(i[0]), 'target': str(i[1]),
                               'color': color, 'penwidth': str(10*math.sqrt(i[2]/max_width))}})

        elif draw_all == True:
            edges_list.append({'data':{'source': str(i[0]), 'target': str(i[1]),
                               'color': 'grey', 'penwidth': str(10*math.sqrt(i[2]/max_width))}})

    # print(edges_list)

    return {'nodes': nodes_list, 'edges': edges_list}



def compute_frequence(edges):
    edges.sort()
    for i in range(len(edges) - 1):
        for j in range(i+1, len(edges)):
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
