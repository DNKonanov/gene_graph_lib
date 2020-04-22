import random
import time
from collections import OrderedDict
import sqlite3
import math
import os

class GenomeGraph:

	"""
	A class used to represent and manipulate a graph form of genomes

	...

	Attributes
	----------
	dict_graph : dict
		a dictionary form of graph representation
	dict_graph_freq : dict
		a dictionary form of graph represenration with frequncy for 
		each edge
	list_graph : dict
		structure which store sequnces of nodes for each contig.
		Access to elements is by genome code and contig code
	genes_code : dict
		table which bind real name of orthology group and its code 
		(real name is key, code is value)
	genes_decode : dict
		table which bind code of orthology group and its real name 
		(code is key, real name is value)

	Methods
	-------
	read_graph(file=None, names_list='all', generate_freq=False)
		Reads a sif file and creates a graph structure

	compute_complexity(outdir, reference, window=20, iterations=500, min_depth=0, max_depth=-1, save_db=None)
		Computes all types of complexity with the full-graph approach

	compute_subgraph_complexity(outdir, reference, window=20, iterations=500, min_depth=0, max_depth=-1, save_db=None)
		Computes all types of complexity with the subgraph approach

	generate_subgraph(self, start_node, end_node, reference, window=20, tails=0, depth=-1, minimal_edge=1)
		Generates subgraph with user-dfined parameters

	find_paths(start, main_chain, min_depth=0, max_depth=-1, window=20)
		Returns deviating paths set generated by "by genomes" method
		
	find_probabilistic_paths(start, main_chain, iterations=500, min_depth=0, max_depth=-1, window=20)
		Returns deviating paths set generated by probabilistic method
	
	find_template(template, method='st', ref='default')
		Find user-defined template in graph structure
	
	find_local_paths(subgr, start, main_chain, window=20)
		Returns deviating paths set in subgraph generated by "by genomes" method

	find_local_probabilistic_paths(subgr, start, main_chain, iterations=500, min_depth=0, max_depth=-1, window=20)
		Returns deviating paths set in subgraph generated by probabilistic method

	save_to_db(self, data, out_db, contig, stamm, window=20)
		Saves complexity data to database

	save_data(self, data, outdir, contig)
		Saves compelxity data in txt format

	"""
	def __init__(self, name='GenomeGraph', file=None):
		"""
		Parameters
		----------
		name : str, optional
			The name of graph structure (default is 'GenomeGraph')
		"""
		self.name = name

	dict_graph = {}
	dict_graph_freq = {}
	list_graph = {}
	genes_code = {}
	genes_decode = {}
	edges_weights = {}
	genes_info = {}

	def _code_genes(self, edge_table):


		all_genes = list(set([edge[0] for edge in edge_table]).union(set([edge[1] for edge in edge_table])))
		self.genes_code = {all_genes[i] : i for i in range(len(all_genes))}
		self.genes_decode = {i : all_genes[i] for i in range(len(all_genes))}


	def _delete_anomaly(self, length_variation=0.05, anomaly_genomes=0.2):

		for node in self.dict_graph:
			if len(self.dict_graph[node]) >= len(self.list_graph)*anomaly_genomes:
				check_anomaly = set([])
				for anode in self.dict_graph[node]:

					if anode not in self.dict_graph:
						continue
					
					if len(self.dict_graph[anode]) > 1:
						break
					check_anomaly.add(self.dict_graph[anode][0])


				if len(check_anomaly) > 1:
					continue

				self.dict_graph[node] = [self.dict_graph[anode][0]]

				for genome in self.list_graph:
					for contig in self.list_graph[genome]:
						c = self.list_graph[genome][contig]

						try:
							if node not in c:
								continue

							pre_index = c.index(node)
							post_index = c.index(self.dict_graph[anode][0])
							if abs(pre_index - post_index):
								c.pop(pre_index + 1)
						except:
							continue



	def read_graph(self, file=None, names_list='all', generate_freq=False):
		"""Reads a sif file and creates a graph structure

		Parameters
		----------
		file : str
			Path to the sif file, generated by orthofinder_parse.py or manually

			Exaple file:
			GENE0001 GENE0002 Genome1 Contig1	(divided by space)
			...
			GENE3999 GENE4000 Genome1 Contig1
			GENE0001 GENE0002 Genome2 Contig2
			...

		names_list : str, optional
			Path to the file with list of interesting genomes

			Example file:
			Genome1
			Genome2
			Genome3
			...

			All genomes from the sif file will be used by default.

	 

		"""
		
		if file == None:
			print('File was not chosen')
			return

		try:
			f_in = open(file, 'r')
		except FileNotFoundError:
			print('File not found!')
			return
		
		edge_table = []

		if names_list != 'all':
			names_list = [name[:-1] for name in open(names_list, 'r')]

		for line in f_in:
			string = line.split(' ')
			string[-1] = string[-1][:-1]

			if len(string) < 3:
				continue

			start_gene = string[0]
			end_gene = string[1]
			stamm = string[2]
			contig = string[3]

			
			if contig not in self.genes_info:
				self.genes_info[contig] = {}

			start_coord = string[4]
			end_coord = string[5]

			self.genes_info[contig][start_gene] = start_coord
			self.genes_info[contig][end_gene] = end_coord
			



			if names_list != 'all' and stamm in names_list:
				edge_table.append([start_gene, end_gene, stamm, contig])
				continue

			elif names_list == 'all':
				edge_table.append([start_gene, end_gene, stamm, contig])

		self._code_genes(edge_table)
		self.list_graph = {edge[2]: {} for edge in edge_table}
		self.dict_graph = {self.genes_code[edge[0]] : set([]) for edge in edge_table}
		if generate_freq != False:
			self.dict_graph_freq = {self.genes_code[edge[0]] : [] for edge in edge_table}


		current_contig = None

		for edge in edge_table:

			start = self.genes_code[edge[0]]
			end = self.genes_code[edge[1]]
			name = edge[2]
			contig = edge[3]
			
			self.dict_graph[start].add(end)
			
			if generate_freq != False:
				self.dict_graph_freq[start].append(end)
			
			if current_contig == None or current_contig != contig:
				self.list_graph[name].update([(contig, [start, end])])
				current_contig = contig

			else:
				self.list_graph[name][contig].append(end)
			
		for gene in self.dict_graph:
			self.dict_graph[gene] = tuple(self.dict_graph[gene])

		if generate_freq != False:

			for gene in self.dict_graph_freq:
				
				freq_table = {}

				for i in self.dict_graph_freq[gene]:
					if i in freq_table:
						continue

					freq_table[i] = self.dict_graph_freq[gene].count(i)

				self.dict_graph_freq[gene] = [(i, freq_table[i]) for i in freq_table]

				for i in freq_table:
					self.edges_weights[(gene, i)] = freq_table[i]

		#self._delete_anomaly()


	def _find_node_env(self, gene, ref):
		
		forward_paths = []
		reversed_paths = []

		for name in self.list_graph:
			for contig in self.list_graph[name]:
				c = self.list_graph[name][contig]
				try:
					index = c.index(gene)
				except: continue
				path = [gene]
				is_uniq = 1
				uniq_counter = 0

				index += 1
				while index < len(c):
					path.append(c[index])
					if is_uniq == 1:
						if c[index] not in ref:
							uniq_counter += 1
						else:
							is_uniq = 0
					index += 1

				forward_paths.append(path + [uniq_counter])

				index = c.index(gene)
				path = [gene]
				is_uniq = 1
				uniq_counter = 0

				index -= 1
				while index >= 0:
					path.append(c[index])
					if is_uniq == 1:
						if c[index] not in ref:
							uniq_counter += 1
						else:
							is_uniq = 0
					index -= 1
					
				reversed_paths.append(path + [uniq_counter])

		return (forward_paths, reversed_paths)

		
		


	def find_template(self, template, method='st', ref='default'):
		"""Find user-defined template in graph structure

		THIS METHOD DOES NOT WORK IN THIS ALPHA VERSION!

		Parameters
		----------
		template : str
			The subgraph template, coded by specific SGL-syntaxis

		method: str, optional
			The search method ('st' - by genomes, 'prob' - probabilistic)
			Default is 'st'

		ref: str, optional
			The reference genome name.
			Default is 'default' and all genomes will be used in the search

		Returns
		-------
		list
			list of possible real nodes names coded by user input 
			['reference_name', {'1': 'OG0000001, '2': 'OG0000055, ...}]
		"""
		print('Input:')
		print(template)
		print()

		template = self._parse_template(template)

		try:
			depth = max([max(path['l']) for path in template['paths']])
		except:
			print('length was chosen automatically')
			depth = 100 

		
		if ref == 'default':
			ref = [r for r in self.list_graph][0]
			contig = [c for c in self.list_graph[ref]][0]


		base = self.list_graph[ref][contig]

		for gene in base:
			print(max(0, base.index(gene) - depth))
			paths = self._find_node_env(gene, base[max(0, base.index(gene) - depth):min(base.index(gene) + depth, len(base))])
			print(len([p for p in paths[0] if p[-1] != 0]))
			print(len([p for p in paths[1] if p[-1] != 0]))

			for edge in template['edges']:
				node_1 = gene

				for node in self.dict_graph_freq[gene]:
					if len(edge['weight']) == 0:
						node_2 = node
					elif node[1] in edge['weigth']:
						node_2 = node
					
		return

		#hits is list of possible real nodes names coded by user input ['reference_name', {'1': 'OG0000001, '2': 'OG0000055, ...}]

	
	def find_paths(self, start, main_chain, min_depth=0, max_depth=-1, window=20):
		"""Returns deviating paths set generated by "by genomes" method

		Parameters
		----------
		start : int
			The name of the node from which the paths will be generated

		main_chain : list or tuple
			The sequence of nodes in the current reference contig

		min_depth : int, optional
			The minimal depth of graph walking
			Default is 0

		max_depth : int, optional
			The maximum depth of graph walking
			Default is unlimited
		
		window : int, optional
			The size of sliding window.
			Default is 20

		Returns
		-------
		list
			List of unique deviating paths generated from start node
			by 'by genomes' method

		"""
		paths = []

		if max_depth == -1:
			max_depth = 100000

		start_index = main_chain.index(start)

		for stamm in self.list_graph:
			for contig in self.list_graph[stamm]:
				stamm_chain = self.list_graph[stamm][contig]

				if start in stamm_chain:
					path = [start]

					index = stamm_chain.index(start)
					for gene in stamm_chain[index + 1:]:
						path.append(gene)
						if len(path) > max_depth:
							break
						if gene in main_chain:

							if abs(start_index - main_chain.index(gene)) > window:
								continue
							if path in paths:
								break
							if len(path) >= min_depth:
								paths.append(path)
							break

		return paths


	def generate_subgraph(self, start_node, end_node, reference, window=20, tails=0, depth=-1, minimal_edge=1):
		"""Generates subgraph with user-dfined parameters

		Because this method is used by GCB server, 
		start_node and end_node are used in decoded form

		Parameters
		----------
		start_node : str
			Real name of the start node in the refernce sequence

		end_node : str
			Real name of the end node in the refernce sequence

		reference : str
			Genome RefSeq name of the reference chain

		window : int, optional
			The number of nodes which will be added to left and right 
			sides of the reference chain
			Default is 20

		depth: int, optional
			The maximum length of stored deviating paths. 
			If path is too long, it will be replaced by left and right tails 
			with length defined in 'tails' parameter
			Default is infinity 
		
		tails : int, optional
			Tails parameter (see in the depth parameter description)
			Default is 0

		minimal_edge : int, optional
			The minimal frequency of edge which will be added to the subgraph
			Default is 1

		Returns
		-------
		dict
			subgraph dictionary representation

		list
			the target chain sequence
			
		"""

		print('Generating subgraph...')
		subgraph = {}
		print('Reference is ' + reference)
		c = 0
		for c in self.list_graph[reference]:
			try:
				start = self.list_graph[reference][c].index(self.genes_code[start_node])
				end = self.list_graph[reference][c].index(self.genes_code[end_node])
				break
				
			except ValueError:
				continue

		if start > end:
			start, end = end, start

		base_chain = self.list_graph[reference][c][max(0, start - window):min(len(self.list_graph[reference][c]), end + window + 1)]

		if depth == -1:
			depth = len(base_chain)

		aim_chain = self.list_graph[reference][c][start: end + 1]

		subgraph[reference] = [base_chain[:]]

		for stamm in self.list_graph:
			if stamm not in subgraph:
				contigs = [[]]

				for contig in self.list_graph[stamm]:
					contigs.append([])
					j = self.list_graph[stamm][contig]
					current_depth = 0
					for gene in range(len(j)):
						if j[gene] in base_chain and contigs[-1] == []:
							current_depth = 0
							contigs[-1] = j[max(0, gene - tails):gene + 1]

						else:
							if contigs[-1] == []:
								continue
							
							contigs[-1].append(j[gene])
							if contigs[-1][-1] in base_chain and contigs[-1][-2] in base_chain:
								current_depth = 0
								continue

							else:
								current_depth += 1
							
							if current_depth >= depth:
								contigs[-1] = contigs[-1][:min(-depth + tails, 0)]
								current_depth = 0
								contigs.append([])
							
							if j[gene] == j[-1]:
								contigs[-1] = contigs[-1][:min(-current_depth + tails, 0)]
								current_depth = 0
								contigs.append([])
								continue
				subgraph[stamm] = contigs

		All_nodes = set([])

		for stamm in subgraph:
			for contig in subgraph[stamm]:
				for gene in contig:
					All_nodes.add(gene)

		i = start-window-1
		while True:
			try:
				if self.list_graph[reference][c][i] in All_nodes:
					subgraph[reference][0] = [self.list_graph[reference][c][i]] + subgraph[reference][0]
					i -= 1
				else:
					break

			except IndexError:
				break

		i = end+window+1
		while True:
			try:
				if self.list_graph[reference][c][i] in All_nodes:
					subgraph[reference][0] = subgraph[reference][0] + [self.list_graph[reference][c][i]]
					i += 1
				else:
					break
			except IndexError:
				break

		print('Completed!')
		

		return subgraph, aim_chain

	
	def find_probabilistic_paths(self, start, main_chain, iterations=500, min_depth=0, max_depth=-1, window=20):
		"""Returns deviating paths set generated by probabilistic method

		Parameters
		----------
		start : int
			The name of the node from which the paths will be generated

		main_chain : list or tuple
			The sequence of nodes in the current reference contig

		iterations : int, optional
			The number of random walks iterations to generate paths

		min_depth : int, optional
			The minimal depth of graph walking
			Default is 0

		max_depth : int, optional
			The maximum depth of graph walking
			Default is unlimited
		
		window : int, optional
			The size of sliding window.
			Default is 20

		Returns
		-------
		list
			List of unique deviating paths generated from start node 
			by probabilistic method

		"""
		paths = []

		if max_depth == -1:
			max_depth = 100000

		start_index = main_chain.index(start)
		for i in range(iterations):
			path = [start]
			current_gene = start
			break_point = 0
			while len(path) <= max_depth:
				if current_gene not in self.dict_graph:
					break

				if break_point > 10:
					break	

				r = random.randint(0, len(self.dict_graph[current_gene]) - 1)
				
				current_gene = self.dict_graph[current_gene][r]

				if current_gene in path:
					break_point += 1
					continue
				break_point = 0
				path.append(current_gene)
								
				if path[-1] in main_chain:
					if path in paths:
						break
					if len(path) >= min_depth:
						paths.append(path)
					break

		return paths


	def _save_to_db(self, data, out_db, contig, stamm, window=20):
		
		db_ = sqlite3.connect(out_db)

		c = db_.cursor()

		contig_id = [q[0] for q in c.execute('select contig_id from contigs_table where contig_code="' + contig + '"')][0]
		nodes = {n[1]: n[0] for n in c.execute('select node_id, node_name from nodes_table where contig_id=' + str(contig_id) + '')}

		for gene in data[0]:
			
			c.execute('insert into complexity_table values(' + str(nodes[self.genes_decode[gene]]) + ', ' + str(contig_id) + ', ' + str(data[0][gene]) + ', ' + 
						str(data[3][gene]) + ', ' + str(data[1][gene]) + ',' + str(data[4][gene]) + ', ' + str(window) + ')')
		db_.commit()
		db_.close()

	def _save_data(self, data, outdir, contig):

		try:
			os.mkdir(outdir + '/extended_info/')
		except:
			pass

		f_io = open(outdir + '/extended_info/IO_complexity_table_contig_' + str(contig) + '.txt', 'a+')
		f_ab = open(outdir + '/extended_info/all_bridges_contig_' + str(contig) + '.txt', 'a+')
		f_wc = open(outdir + '/window_complexity_contig_' + str(contig) + '.txt', 'a+')
		f_mc = open(outdir + '/extended_info/main_chain_contig_' + str(contig) + '.txt', 'a+')

		f_wc.write('position\tOrthoGroupID\tcomplexity\n')

		for gene in data[0]:
			f_wc.write(self.genes_info[contig][self.genes_decode[gene]] + '\t' + self.genes_decode[gene] + '\t' + str(data[0][gene]) + '\n')
			f_mc.write(self.genes_info[contig][self.genes_decode[gene]] + '\t' + self.genes_decode[gene] + '\n')

		for gene in data[1]:
			f_io.write(self.genes_info[contig][self.genes_decode[gene]] + '\t' + self.genes_decode[gene] + '\t' + str(data[1][gene]) + '\n')

		for bridge in data[2]:
			f_ab.write(self.genes_decode[bridge[0]] + '\t' + self.genes_decode[bridge[1]] + '\t' + str(data[2][bridge]) + '\n')

		f_io.close()
		f_ab.close()
		f_wc.close()
		f_mc.close()
		
		f_io = open(outdir + '/extended_info/prob_IO_complexity_table_contig_' + str(contig) + '.txt', 'a+')
		f_ab = open(outdir + '/extended_info/prob_all_bridges_contig_' + str(contig) + '.txt', 'a+')
		f_wc = open(outdir + '/prob_window_complexity_contig_' + str(contig) + '.txt', 'a+')
		f_mc = open(outdir + '/extended_info/prob_main_chain_contig_' + str(contig) + '.txt', 'a+')
		
		f_wc.write('position\tOrthoGroupID\tcomplexity\n')

		for gene in data[3]:
			f_wc.write(self.genes_info[contig][self.genes_decode[gene]] + '\t' + self.genes_decode[gene] + '\t' + str(data[3][gene]) + '\n')
			f_mc.write(self.genes_info[contig][self.genes_decode[gene]] + '\t' + self.genes_decode[gene] + '\n')

		for gene in data[4]:
			f_io.write(self.genes_info[contig][self.genes_decode[gene]] + '\t' + self.genes_decode[gene] + '\t' + str(data[4][gene]) + '\n')

		for bridge in data[5]:
			f_ab.write(self.genes_decode[bridge[0]] + '\t' + self.genes_decode[bridge[1]] + '\t' + str(data[5][bridge]) + '\n')

		f_io.close()
		f_ab.close()
		f_wc.close()
		f_mc.close()

		
	def compute_complexity(self, outdir, reference, window=20, iterations=500, min_depth=0, max_depth=-1, save_db=None):
		"""Computes all types of complexity with the full-graph approach

		Parameters
		----------
		outdir : str
			The name of output directory.

		reference : str
			The name of the reference genome

		window : int, optional
			The size of sliding window
			Default is 20

		iterations : int, optional
			The number of random walks iterations to generate paths
			in propabilistic method
			Default is 500

		min_depth : int, optional
			The minimal depth of graph walking
			Default is 0

		max_depth : int, optional
			The maximum depth of graph walking
			Default is unlimited

		save_db : bool or str
			If save_d parameter is not None, all complexity data will be saved
			to user-defined database

		"""
		
		print('Reference is ' + reference)
		print('Number of contigs: ' + str(len(self.list_graph[reference])))


		for contig in self.list_graph[reference]:
			print('\nComputing wth contig ' + str(contig) + '...')

			base_line = self.list_graph[reference][contig]
			
			complexity_table = OrderedDict((gene, 0) for gene in base_line)
			prob_complexity_table = OrderedDict((gene, 0) for gene in base_line)
			io_table = OrderedDict((gene, 0) for gene in base_line)
			prob_io_table = OrderedDict((gene, 0) for gene in base_line)
			all_bridges = OrderedDict([])
			prob_all_bridges = OrderedDict([])
			norm = 2*window
			count = 1
			for gene in base_line:
				print(str(count) + ' gene of ' + str(len(base_line)), end='')

				
				#by stamm method
				paths = self.find_paths(gene, base_line, min_depth=min_depth, max_depth=max_depth, window=window)
				for p in paths:

					if (p[0], p[-1]) not in all_bridges:

						all_bridges[(p[0], p[-1])] = 1
					else:
						all_bridges[(p[0], p[-1])] += 1
					
					io_table[p[0]] += 1
					io_table[p[-1]] += 1

					start_index, end_index = min(base_line.index(p[0]), base_line.index(p[-1])), max(base_line.index(p[0]), base_line.index(p[-1]))
					
					if abs(start_index - end_index) <= window:
						if len(p) == 2 and end_index-start_index == 1:
							continue
						for i in range(start_index, end_index + 1):
							complexity_table[base_line[i]] += 1/norm
	
				
				
				#probabilistic method
				paths = self.find_probabilistic_paths(gene, base_line, iterations=iterations, min_depth=min_depth, max_depth=max_depth, window=window)

				for p in paths:

					if (p[0], p[-1]) not in prob_all_bridges:

						prob_all_bridges[(p[0], p[-1])] = 1
					else:
						prob_all_bridges[(p[0], p[-1])] += 1
					
					prob_io_table[p[0]] += 1
					prob_io_table[p[-1]] += 1
					
					start_index, end_index = min(base_line.index(p[0]), base_line.index(p[-1])), max(base_line.index(p[0]), base_line.index(p[-1]))

					if abs(start_index - end_index) <= window:
						if len(p) == 2 and end_index-start_index == 1:
							continue
						for i in range(start_index, end_index + 1):
							prob_complexity_table[base_line[i]] += 1/norm

				print('\r', end='')
				count += 1

			self._save_data([complexity_table, io_table, all_bridges, 
						prob_complexity_table, prob_io_table, prob_all_bridges], outdir, contig)

			if save_db != None:
				self._save_to_db([complexity_table, io_table, all_bridges, 
						prob_complexity_table, prob_io_table, prob_all_bridges], save_db, contig, reference, window=window)
	
		print('\nComputing completed')



	def generate_local_subgraph(self, start_node, end_node, reference, window=20, depth=-1, minimal_edge=1):
		"""Generates subgraph to compute subgraph complexity

		Because this method is used by GCB server, 
		start_node and end_node are used in decoded form

		Parameters
		----------
		start_node : str
			Real name of the start node in the refernce sequence

		end_node : str
			Real name of the end node in the refernce sequence

		reference : str
			Genome RefSeq name of the reference chain

		window : int, optional
			The number of nodes which will be added to left and right 
			sides of the reference chain
			Default is 20

		depth: int, optional
			The maximum length of stored deviating paths. 
			If path is too long, it will be replaced by left and right tails 
			with length defined in 'tails' parameter
			Default is infinity 
		
		minimal_edge : int, optional
			The minimal frequency of edge which will be added to the subgraph
			Default is 1

		Returns
		-------
		dict
			subgraph dictionary representation

		dict
			subgraph list representation
			
		"""
		subgraph = {}
		c = 0
		for c in self.list_graph[reference]:
			try:
				start = self.list_graph[reference][c].index(start_node)
				end = self.list_graph[reference][c].index(end_node)
				break
				
			except ValueError:
				continue

		if start > end:
			start, end = end, start

		base_chain = self.list_graph[reference][c][max(0, start - window):min(len(self.list_graph[reference][c]), end + window + 1)]

		if depth == -1:
			depth = len(base_chain)

		aim_chain = self.list_graph[reference][c][start: end + 1]

		subgraph[reference] = [base_chain[:]]

		for stamm in self.list_graph:
			if stamm not in subgraph:
				contigs = [[]]

				for contig in self.list_graph[stamm]:
					contigs.append([])
					j = self.list_graph[stamm][contig]
					current_depth = 0
					for gene in range(len(j)):
						if j[gene] in base_chain and contigs[-1] == []:
							current_depth = 0
							contigs[-1] = j[max(0, gene):gene + 1]

						else:
							if contigs[-1] == []:
								continue
							
							contigs[-1].append(j[gene])
							if contigs[-1][-1] in base_chain and contigs[-1][-2] in base_chain:
								current_depth = 0
								continue

							else:
								current_depth += 1
							
							if current_depth >= depth:
								contigs[-1] = contigs[-1][:min(-depth, 0)]
								current_depth = 0
								contigs.append([])
							
							if j[gene] == j[-1]:
								contigs[-1] = contigs[-1][:min(-current_depth, 0)]
								current_depth = 0
								contigs.append([])
								continue
				subgraph[stamm] = contigs

		All_nodes = set([])


		for stamm in subgraph:
			for contig in subgraph[stamm]:
				for gene in contig:
					All_nodes.add(gene)

		i = start-window-1
		while True:
			try:
				if self.list_graph[reference][c][i] in All_nodes:
					subgraph[reference][0] = [self.list_graph[reference][c][i]] + subgraph[reference][0]
					i -= 1
				else:
					break

			except IndexError:
				break

		i = end+window+1
		while True:
			try:
				if self.list_graph[reference][c][i] in All_nodes:
					subgraph[reference][0] = subgraph[reference][0] + [self.list_graph[reference][c][i]]
					i += 1
				else:
					break
			except IndexError:
				break
		dict_graph = {}
		for name in subgraph:
			for contig in subgraph[name]:
				for i in range(len(contig) - 1):
					if contig[i] not in dict_graph:
						dict_graph[contig[i]] = set([contig[i+1]])
					else:
						dict_graph[contig[i]].add(contig[i+1])


		for gene in dict_graph:
			dict_graph[gene] = tuple(dict_graph[gene])
		
		return dict_graph, subgraph

	def find_local_paths(self, subgr, start, main_chain, window=20):
		"""Returns deviating paths set generated by probabilistic method.

		This method is used to compute subgraph complexity

		Parameters
		----------
		subgr : dict
			List form of the input subgraph.
			Generated by generate_local_subgraph method

		start : int
			The name of the node from which the paths will be generated

		main_chain : list or tuple
			The sequence of nodes in the current reference contig

		window : int, optional
			The size of sliding window.
			Default is 20

		Returns
		-------
		list
			List of unique deviating paths generated from start node 
			by 'by genomes' method in the input subgraph

		"""
		paths = []

		start_index = main_chain.index(start)
		for stamm in subgr:
			for contig in subgr[stamm]:
				stamm_chain = contig

				if start in stamm_chain:
					path = [start]

					index = stamm_chain.index(start)
					for gene in stamm_chain[index + 1:]:
						path.append(gene)
						if gene in main_chain:

							if abs(start_index - main_chain.index(gene)) > window:
								continue
							if path in paths:
								break
							paths.append(tuple(path))
							break

		return paths


	def find_local_probabilistic_paths(self, subgr, start, main_chain, iterations=500, min_depth=0, max_depth=-1, window=20):
		"""Returns deviating paths set generated by probabilistic method.

		This method is used to compute subgraph complexity

		Parameters
		----------
		subgr : dict
			Dictionary form of the input subgraph.
			Generated by generate_local_subgraph method

		start : int
			The name of the node from which the paths will be generated

		main_chain : list or tuple
			The sequence of nodes in the current reference contig

		iterations : int, optional
			The number of random walks iterations to generate paths

		min_depth : int, optional
			The minimal depth of graph walking
			Default is 0

		max_depth : int, optional
			The maximum depth of graph walking
			Default is unlimited
		
		window : int, optional
			The size of sliding window.
			Default is 20

		Returns
		-------
		list
			List of unique deviating paths generated from start node 
			by probabilistic method in the input subgraph

		"""
		paths = []

		if max_depth == -1:
			max_depth = 100000

		start_index = main_chain.index(start)
		for i in range(iterations):
			path = [start]
			current_gene = start
			break_point = 0
			while len(path) <= max_depth:
				if current_gene not in subgr:
					break

				if break_point > 10:
					break	

				r = random.randint(0, len(subgr[current_gene]) - 1)
				
				current_gene = subgr[current_gene][r]

				if current_gene in path:
					break_point += 1
					continue
				break_point = 0
				path.append(current_gene)
								
				if path[-1] in main_chain:
					if path in paths:
						break
					if len(path) >= min_depth:
						paths.append(tuple(path))
					break
		return paths

			
	def compute_subgraph_complexity(self, outdir, reference, window=20, iterations=500, min_depth=0, max_depth=-1, save_db=None):
		"""Computes all types of complexity with the full-graph approach

		Parameters
		----------
		outdir : str
			The name of output directory.

		reference : str
			The name of the reference genome

		window : int, optional
			The size of sliding window
			Default is 20

		iterations : int, optional
			The number of random walks iterations to generate paths
			in propabilistic method
			Default is 500

		min_depth : int, optional
			The minimal depth of graph walking
			Default is 0

		max_depth : int, optional
			The maximum depth of graph walking
			Default is unlimited

		save_db : bool or str
			If save_d parameter is not None, all complexity data will be saved
			to user-defined database

		"""
		
		print('Reference is ' + reference)
		print('Number of contigs: ' + str(len(self.list_graph[reference])))


		for contig in self.list_graph[reference]:
			print('\nComputing wth contig ' + str(contig) + '...')

			base_line = self.list_graph[reference][contig]
			
			complexity_table = OrderedDict((gene, 0) for gene in base_line)
			prob_complexity_table = OrderedDict((gene, 0) for gene in base_line)
			io_table = OrderedDict((gene, 0) for gene in base_line)
			prob_io_table = OrderedDict((gene, 0) for gene in base_line)
			all_bridges = OrderedDict([])
			prob_all_bridges = OrderedDict([])
			norm = 2*window
			count = 1

			All_paths = set([])
			for gene in base_line:
				print(str(count) + ' gene of ' + str(len(base_line)), end='')

				
				subgr, list_subgr = self.generate_local_subgraph(gene, gene, reference, window=window, depth=max_depth)


				for local_gene in list_subgr[reference][0]:
					

					#by genomes method with subgraph
					paths = self.find_local_paths(list_subgr, local_gene, list_subgr[reference][0], window=window)

					for p in paths:
						if p in All_paths:
							continue

						if (p[0], p[-1]) not in all_bridges:

							all_bridges[(p[0], p[-1])] = 1
						else:
							all_bridges[(p[0], p[-1])] += 1
						
						io_table[p[0]] += 1
						io_table[p[-1]] += 1
						
						start_index, end_index = min(base_line.index(p[0]), base_line.index(p[-1])), max(base_line.index(p[0]), base_line.index(p[-1]))
						All_paths.add(p)

						if abs(start_index - end_index) <= window:
							if len(p) == 2 and end_index-start_index == 1:
								continue
							for i in range(start_index, end_index + 1):
								complexity_table[base_line[i]] += 1/norm

					#probabilistic method with subgraph
					paths = self.find_local_probabilistic_paths(subgr, local_gene, list_subgr[reference][0], iterations=iterations, min_depth=min_depth, max_depth=max_depth, window=window)

					for p in paths:
						if p in All_paths:
							continue

						if (p[0], p[-1]) not in prob_all_bridges:

							prob_all_bridges[(p[0], p[-1])] = 1
						else:
							prob_all_bridges[(p[0], p[-1])] += 1
						
						prob_io_table[p[0]] += 1
						prob_io_table[p[-1]] += 1
						
						start_index, end_index = min(base_line.index(p[0]), base_line.index(p[-1])), max(base_line.index(p[0]), base_line.index(p[-1]))
						All_paths.add(p)

						if abs(start_index - end_index) <= window:
							if len(p) == 2 and end_index-start_index == 1:
								continue
							for i in range(start_index, end_index + 1):
								prob_complexity_table[base_line[i]] += 1/norm

				print('\r', end='')
				count += 1

			self._save_data([complexity_table, io_table, all_bridges, 
						prob_complexity_table, prob_io_table, prob_all_bridges], outdir, contig)

			if save_db != None:
				self._save_to_db([complexity_table, io_table, all_bridges, 
						prob_complexity_table, prob_io_table, prob_all_bridges], save_db, contig, reference, window=window)
	
		print('\nComputing completed')

