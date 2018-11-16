import random
import time
from collections import OrderedDict
import sqlite3
import math

class GenomeGraph:
	def __init__(self, name='GenomeGraph', file=None):
		...
		self.name = name

	dict_graph = {}
	dict_graph_freq = {}
	list_graph = {}
	genes_code = {}
	genes_decode = {}


	def code_genes(self, edge_table):

		all_genes = list(set([edge[0] for edge in edge_table]).union(set([edge[1] for edge in edge_table])))
		self.genes_code = {all_genes[i] : i for i in range(len(all_genes))}
		self.genes_decode = {i : all_genes[i] for i in range(len(all_genes))}




	def read_graph(self, file=None, names_list='all', generate_freq=False):
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

			if len(string) < 3:
				continue

			start_gene = string[0]
			end_gene = string[1]
			stamm = string[2]
			contig = string[3][:-1]

			if names_list != 'all' and stamm in names_list:
				edge_table.append([start_gene, end_gene, stamm, contig])
				continue

			elif names_list == 'all':
				edge_table.append([start_gene, end_gene, stamm, contig])

		self.code_genes(edge_table)
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


	def _parse_template(self, template):

		parsed_template = {'nodes': [], 'edges': [], 'paths': []}

		skip_sb = [' ', '\n', '\t']
		
		for sb in skip_sb:
			template.replace(sb, '')

		open_node = 'close'
		i = 0
		while i < len(template):
			sb = template[i]
			if sb == "'":
				i += 1
				string = ''
				while template[i] != "'":
					string += template[i]
					i += 1
				
				node = self._parse_node(string)
				parsed_template['nodes'].append(node)
				if open_node == 'edge':
					parsed_template['edges'][-1]['end'] = node
					open_node = 'close'
				if open_node == 'path':
					parsed_template['paths'][-1]['end'] = node
					open_node = 'close'
				
				i += 1
				continue

			if sb == '-':
				i += 1
				string = ''
				while template[i] != '>':
					string += template[i]
					i += 1

				try:
					edge = self._parse_edge(string, parsed_template['nodes'][-1])

				except IndexError:
					edge = self._parse_edge(string, '')

				parsed_template['edges'].append(edge)
				open_node = 'edge'
				i += 1
				continue

			if sb == '*':
				i += 1
				string = ''
				while template[i] != "/":
					string += template[i]
					i += 1

				try:
					path = self._parse_path(string, parsed_template['nodes'][-1])
				except IndexError:
					path = self._parse_path(string, '')
				parsed_template['paths'].append(path)
				open_node = 'path'

			if sb == ';':
				open_node = 'close'
				i += 1
				continue
			i += 1

		parsed_template['nodes'] = set(parsed_template['nodes'])

		return parsed_template

	def _parse_node(self, string):
		node = string
		return node

	def _parse_edge(self, string, last_node):
		start = last_node
		weigth = set([])
		end = ''


		i = 0
		if string.startswith('{'):
			string = string[1:-1]
			parameters = string.split(':')
			for par in parameters:
				if par.startswith('!') == False:
					if par.startswith('w('):
						i, j = int(par[2:-1].split(',')[0]), int(par[2:-1].split(',')[1])
						weigth = weigth.union(set(list(weigth) + [n for  n in range(i, j+1)]))
			for par in parameters:
				if par.startswith('!'):
					if par.startswith('w('):
						i, j = int(par[3:-1].split(',')[0]), int(par[2:-1].split(',')[1])
						weigth = weigth.difference(set([n for  n in range(i, j+1)]))

		string = {'start': start, 'end': end, 'w': weigth}

		return string


	def _parse_path(self, string, last_node):
		start = last_node
		weigth = set([])
		uniq = set([])
		n = []
		nr = False
		length = set([])
		weigth
		end = ''

		i = 0
		if string.startswith('{'):
			string = string[1:-1]
			parameters = string.split(':')
			for par in parameters:
				if par.startswith('!') == False:
					if par.startswith('l('):
						i, j = int(par[2:-1].split(',')[0]), int(par[2:-1].split(',')[1])
						length = length.union(set(list(length) + [n for  n in range(i, j+1)]))
					if par.startswith('w('):
						i, j = int(par[2:-1].split(',')[0]), int(par[2:-1].split(',')[1])
						weigth = weigth.union(set(list(weigth) + [n for  n in range(i, j+1)]))
					if par.startswith('u('):
						i, j = int(par[2:-1].split(',')[0]), int(par[2:-1].split(',')[1])
						uniq = uniq.union(set(list(uniq) + [n for  n in range(i, j+1)]))

			
			
			for par in parameters:
				if par.startswith('!'):
					if par.startswith('!l('):
						i, j = int(par[3:-1].split(',')[0]), int(par[2:-1].split(',')[1])
						length = length.difference(set([n for  n in range(i, j+1)]))
					if par.startswith('!w('):
						i, j = int(par[3:-1].split(',')[0]), int(par[2:-1].split(',')[1])
						weigth = weigth.difference(set([n for  n in range(i, j+1)]))
					if par.startswith('!u('):
						i, j = int(par[3:-1].split(',')[0]), int(par[2:-1].split(',')[1])
						uniq = uniq.difference(set([n for  n in range(i, j+1)]))

		string = {'start': start, 'end': end, 'l': length, 'w': weigth, 'u': uniq}

		return string

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

	
	def find_paths(self, start, main_chain, min_depth=0, max_depth=-1):
		paths = []

		if max_depth == -1:
			max_depth = 100000

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
							if path in paths:
								break
							if len(path) >= min_depth:
								paths.append(path)
							break

		return paths


	def generate_subgraph(self, start_node, end_node, reference, window=20, tails=0, depth=-1, minimal_edge=1):

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
								contigs
							
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





	def find_probabilistic_paths(self, start, main_chain, iterations=500, min_depth=0, max_depth=-1):
		paths = []

		if max_depth == -1:
			max_depth = 100000

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


	def save_to_db(self, data, out_db, contig, stamm, window=20):



		connect = sqlite3.connect(out_db)

		c = connect.cursor()


		contig_key = contig_key = [row for row in c.execute('SELECT id FROM contigs_table WHERE  contig = "' + contig + '"')][0][0]


		sql_exec = 'CREATE TABLE IF NOT EXISTS og_complexity_table_' + str(window) + ' (contig integer,og text,win_var float,prob_win_var float,io float,prob_io float	);'
		for gene in data[0]:
			sql_exec = 'INSERT INTO og_complexity_table_' + str(window) + ' VALUES (' + str(contig_key) + ',"' + self.genes_decode[gene] + '",' + str(data[0][gene]) + ',' + str(data[3][gene]) + ',' + str(data[1][gene]) + ',' + str(data[4][gene]) + ')'
			c.execute(sql_exec)

		connect.commit()
		connect.close()
		


	def save_data(self, data, outdir, contig):
		f_io = open(outdir + '/IO_vaiability_table_contig_' + str(contig) + '.txt', 'a+')
		f_ab = open(outdir + '/all_bridges_contig_' + str(contig) + '.txt', 'a+')
		f_wc = open(outdir + '/window_variability_contig_' + str(contig) + '.txt', 'a+')
		f_mc = open(outdir + '/main_chain_contig_' + str(contig) + '.txt', 'a+')

		for gene in data[0]:
			f_wc.write(self.genes_decode[gene] + '\t' + str(data[0][gene]) + '\n')
			f_mc.write(self.genes_decode[gene] + '\n')

		for gene in data[1]:
			f_io.write(self.genes_decode[gene] + '\t' + str(data[1][gene]) + '\n')

		for bridge in data[2]:
			f_ab.write(self.genes_decode[bridge[0]] + '\t' + self.genes_decode[bridge[1]] + '\t' + str(data[2][bridge]) + '\t')

		
		f_io = open(outdir + '/prob_IO_vaiability_table_contig_' + str(contig) + '.txt', 'a+')
		f_ab = open(outdir + '/prob_all_bridges_contig_' + str(contig) + '.txt', 'a+')
		f_wc = open(outdir + '/prob_window_variability_contig_' + str(contig) + '.txt', 'a+')
		f_mc = open(outdir + '/prob_main_chain_contig_' + str(contig) + '.txt', 'a+')

		for gene in data[3]:
			f_wc.write(self.genes_decode[gene] + '\t' + str(data[3][gene]) + '\n')
			f_mc.write(self.genes_decode[gene] + '\n')

		for gene in data[4]:
			f_io.write(self.genes_decode[gene] + '\t' + str(data[4][gene]) + '\n')

		for bridge in data[5]:
			f_ab.write(self.genes_decode[bridge[0]] + '\t' + self.genes_decode[bridge[1]] + '\t' + str(data[5][bridge]) + '\t')


	def compute_variability(self, outdir, reference, window=20, iterations=500, min_depth=0, max_depth=-1, save_db=None):
		
		print('Reference is ' + reference)
		print('Number of contigs: ' + str(len(self.list_graph[reference])))


		for contig in self.list_graph[reference]:
			print('\nComputing wth contig ' + str(contig) + '...')

			base_line = self.list_graph[reference][contig]
			
			variability_table = OrderedDict((gene, 0) for gene in base_line)
			prob_variability_table = OrderedDict((gene, 0) for gene in base_line)
			io_table = OrderedDict((gene, 0) for gene in base_line)
			prob_io_table = OrderedDict((gene, 0) for gene in base_line)
			all_bridges = OrderedDict([])
			prob_all_bridges = OrderedDict([])

			count = 1
			for gene in base_line:
				print(str(count) + ' gene of ' + str(len(base_line)), end='')

				norm = min(len(base_line), base_line.index(gene) + window) - max(0, base_line.index(gene) - window)
				#by stamm method
				paths = self.find_paths(gene, base_line, min_depth=min_depth, max_depth=max_depth)
				for p in paths:

					if (p[0], p[-1]) not in all_bridges:

						all_bridges[p[0], p[-1]] = 1
					else:
						all_bridges[p[0], p[-1]] += 1
					
					io_table[p[0]] += 1
					io_table[p[-1]] += 1
					
					start_index = base_line.index(p[0])
					end_index = base_line.index(p[-1])

					if abs(start_index - end_index) <= window:
						for i in range(min(start_index, end_index), max(start_index, end_index) + 1):
							variability_table[base_line[i]] += 1/float(norm)
				
				
				
				#probabilistic method
				paths = self.find_probabilistic_paths(gene, base_line, iterations=iterations, min_depth=min_depth, max_depth=max_depth)

				for p in paths:

					if (p[0], p[-1]) not in prob_all_bridges:

						prob_all_bridges[p[0], p[-1]] = 1
					else:
						prob_all_bridges[p[0], p[-1]] += 1
					
					prob_io_table[p[0]] += 1
					prob_io_table[p[-1]] += 1

					start_index = base_line.index(p[0])
					end_index = base_line.index(p[-1])

					if abs(start_index - end_index) <= window:
						for i in range(min(start_index, end_index), max(start_index, end_index) + 1):
							prob_variability_table[base_line[i]] += 1/float(norm)

				print('\r', end='')
				count += 1

			self.save_data([variability_table, io_table, all_bridges, 
						prob_variability_table, prob_io_table, prob_all_bridges], outdir, contig)

			if save_db != None:
				self.save_to_db([variability_table, io_table, all_bridges, 
						prob_variability_table, prob_io_table, prob_all_bridges], save_db, contig, reference, window=window)
	
		print('\nComputing completed')
		