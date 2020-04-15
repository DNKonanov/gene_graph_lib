

def find_context(graph):

    og_set = set([])
    edge_og = set([])
    for name in graph:
        for contig in graph[name]:
            c = graph[name][contig]
            

            edge_og.add(c[0])
            edge_og.add(c[-1])
            for og in c:
                og_set.add(og)

    og_set.difference_update(edge_og)

    og_context = {og: [] for og in og_set}

    for name in graph:
        for contig in graph[name]:
            c = graph[name][contig]

            for i in range(1, len(c) - 1):
                if c[i] not in og_context:
                    continue

                if ((c[i-1], c[i+1]) in og_context[c[i]] or (c[i+1], c[i-1]) in og_context[c[i]]):
                    pass
                else:
                    og_context[c[i]].append((c[i-1], c[i+1]))

    return {og: len(og_context[og]) for og in og_context}

    