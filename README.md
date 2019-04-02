# gene_graph_lib

This is a Python module containing methods used in the [GCB](gcb.rcpcm.org) web-application.
Contains classes and methods used to create and process the graph representation of a genomes set. 

## compute_complexity module

This module describes the main object `GenomeGraph`.

### Attributes

attribute | type | description
--------- | ---- | -----------
`dict_graph` | `dict` | a dictionary form of graph representation
`dict_graph_freq` | `dict` | a dictionary form of graph represenration with frequncy for each edge
`list_graph` | `dict` | structure which store sequnces of nodes for each contig. Access to elements is by genome code and contig code
`genes_code` | `dict` | table which bind real name of orthology group and its code (real name is key, code is value)
`genes_decode` | `dict` | table which bind code of orthology group and its real name (code is key, real name is value)


### Methods

method | description | return
------ | ----------- | ------
`read_graph (file=None, names_list='all', generate_freq=False)` | Reads a sif file and creates a graph structure | no returned value
`compute_complexity (outdir, reference, window=20, iterations=500, min_depth=0, max_depth=-1, save_db=None)` | Computes all types of complexity with the full-graph approach | no returned value
`compute_subgraph_complexity (outdir, reference, window=20, iterations=500, min_depth=0, max_depth=-1, save_db=None)` | Computes all types of complexity with the subgraph approach | no returned value
`generate_subgraph (self, start_node, end_node, reference, window=20, tails=0, depth=-1, minimal_edge=1)` | Generates subgraph with user-dfined parameters | subgraph
`find_paths (start, main_chain, min_depth=0, max_depth=-1, window=20)` | Returns deviating paths set generated by "by genomes" method | list of paths
`find_probabilistic_paths (start, main_chain, iterations=500, min_depth=0, max_depth=-1, window=20)` | Returns deviating paths set generated by probabilistic method | list of paths
`find_template (template, method='st', ref='default')` | Find user-defined template in graph structure | template hits list
`find_local_paths (subgr, start, main_chain, window=20)` | Returns deviating paths set in subgraph generated by "by genomes" method | list of paths
`find_local_probabilistic_paths (subgr, start, main_chain, iterations=500, min_depth=0, max_depth=-1, window=20)` | Returns deviating paths set in subgraph generated by probabilistic method | list of paths
`save_to_db (self, data, out_db, contig, stamm, window=20)` | Saves complexity data to database | no returned value
`save_data (self, data, outdir, contig)` | Saves compelxity data in txt format | no returned value

