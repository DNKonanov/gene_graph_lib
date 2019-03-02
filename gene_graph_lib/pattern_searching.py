import networkx
import collections
import compute_complexity
import time


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


def count_garlic_in_gene(genome_dict, garlic_pattern=GarlicPattern()) -> dict:
    result = {}
    hashed = set()
    for genome, contigs in genome_dict.items():
        result[genome] = []
        for contig_name in contigs:
            current_contig = contigs[contig_name]
            for other_genome in genome_dict.keys() - set(genome):
                for other_contig in genome_dict[other_genome].values():
                    result[genome].extend(garlic_number(
                        current_contig,
                        other_contig,
                        min_insert_per_pos=garlic_pattern.min_insert_per_pos,
                        max_insert_per_pos=garlic_pattern.max_insert_per_pos,
                        hashed=hashed))
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

test_graph = compute_complexity.GenomeGraph()
test_graph.read_graph('../Streptococcus_pneumoniae.sif')
print('readed')
# for k, v in test_graph.list_graph.items():
#     print(k, list(map(len, v.values())))

keys_ = [k for k in test_graph.list_graph.keys()][:15]

t = time.time()
for k, v in count_garlic_in_gene({k: test_graph.list_graph[k] for k in keys_}).items():
    print(k, len(v))
print(time.time() - t)
