import random
import time
from collections import OrderedDict


#LEGACY
def divide_by_stamm(edge_table):

	print('\n\tDividing to stamms...')
	stamm_table = OrderedDict([])

	for edge in edge_table:

		start_gene = edge[0]
		end_gene = edge[1]
		stamm = edge[2]

		if stamm not in stamm_table:
			stamm_table.update([(stamm, [[start_gene, end_gene]])])
			continue

		if start_gene == stamm_table[stamm][-1][-1]:
			stamm_table[stamm][-1].append(end_gene)
			continue

		else:
			stamm_table[stamm].append([start_gene, end_gene])
			continue


	print('\tDividing to stamms completed!!!')
	
	return stamm_table


#LEGACY
def generate_graph(edge_table):

	print('-----------------------')
	t = time.time()
	print('Generating of graph...')
	graph = divide_by_stamm(edge_table)

	for stamm in graph:
		graph[stamm].sort(key=len, reverse=True)

	print('Completed!')
	print('Graph generation time:\t' + str(int(time.time()*100 - t*100)/100) + 's\n')

	return graph





class GenomeGraph:
	def __init__(self, name='GenomeGraph', file=None):
		...
		self.name = name

	dict_graph = {}
	list_graph = {}
	genes_code = {}
	genes_decode = {}


	def code_genes(self, edge_table):

		all_genes = list(set([edge[0] for edge in edge_table]).union(set([edge[1] for edge in edge_table])))
		self.genes_code = {all_genes[i] : i for i in range(len(all_genes))}
		self.genes_decode = {i : all_genes[i] for i in range(len(all_genes))}


	def read_graph(self, file=None, names_list='all'):
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
			stamm = string[2][:-1]

			if names_list != 'all' and stamm in names_list:
				edge_table.append([start_gene, end_gene, stamm])
				continue

			elif names_list == 'all':
				edge_table.append([start_gene, end_gene, stamm])

		self.code_genes(edge_table)
		self.list_graph = {edge[2]: [] for edge in edge_table}
		self.dict_graph = {self.genes_code[edge[0]] : set([]) for edge in edge_table}

		for edge in edge_table:

			start = self.genes_code[edge[0]]
			end = self.genes_code[edge[1]]
			name = edge[2]

			self.dict_graph[start].add(end)

			if self.list_graph[name] == []:
				self.list_graph[name].append([start, end])

			elif self.list_graph[name][-1][-1] != start:
				self.list_graph[name].append([start, end])

			else:
				self.list_graph[name][-1].append(end)
			
		for gene in self.dict_graph:
			self.dict_graph[gene] = tuple(self.dict_graph[gene])

	
	def find_paths(self, start, main_chain, min_depth=0, max_depth=-1):
		paths = []

		if max_depth == -1:
			max_depth = 100000

		for stamm in self.list_graph:
			for contig in self.list_graph[stamm]:
				stamm_chain = contig

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


	def generate_subgraph(self, start_node, end_node, reference, window=20, tails=0, depth=-1):

		print('Generating subgraph...')
		subgraph = {}
		print('Reference is ' + reference)
		c = 0
		while(True):
			try:
				start = self.list_graph[reference][c].index(self.genes_code[start_node])
				end = self.list_graph[reference][c].index(self.genes_code[end_node])
				break
				
			except ValueError:
				c += 1
				pass

		if start > end:
			start, end = end, start

		base_chain = self.list_graph[reference][c][max(0, start - window):min(len(self.list_graph[reference][c]) - 1, end + window + 1)]

		if depth == -1:
			depth = len(base_chain)

		aim_chain = self.list_graph[reference][c][start: end + 1]

		subgraph[reference] = [base_chain[:]]

		for stamm in self.list_graph:
			if stamm not in subgraph:
				contigs = [[]]

				for j in self.list_graph[stamm]:
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
								contigs[-1] = contigs[-1][:min(-current_depth + tails, -0)]
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
			if self.list_graph[reference][c][i] in All_nodes:
				subgraph[reference][0] = [self.list_graph[reference][c][i]] + subgraph[reference][0]
				i -= 1
			else:
				break

		i = end+window+1
		while True:
			if self.list_graph[reference][c][i] in All_nodes:
				subgraph[reference][0] = subgraph[reference][0] + [self.list_graph[reference][c][i]]
				i += 1
			else:
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


	def save_data(self, data, outdir, contig):
		f_io = open(outdir + '/IO_vaiability_table_contig' + str(contig) + '.txt', 'a+')
		f_ab = open(outdir + '/all_bridges_contig' + str(contig) + '.txt', 'a+')
		f_wc = open(outdir + '/window_variability_contig' + str(contig) + '.txt', 'a+')
		f_mc = open(outdir + '/main_chain_contig' + str(contig) + '.txt', 'a+')

		for gene in data[0]:
			f_wc.write(self.genes_decode[gene] + '\t' + str(data[0][gene]) + '\n')
			f_mc.write(self.genes_decode[gene] + '\n')

		for gene in data[1]:
			f_io.write(self.genes_decode[gene] + '\t' + str(data[1][gene]) + '\n')

		for bridge in data[2]:
			f_ab.write(self.genes_decode[bridge[0]] + '\t' + self.genes_decode[bridge[1]] + '\t' + str(data[2][bridge]) + '\t')

		
		f_io = open(outdir + '/prob_IO_vaiability_table_contig' + str(contig) + '.txt', 'a+')
		f_ab = open(outdir + '/prob_all_bridges_contig' + str(contig) + '.txt', 'a+')
		f_wc = open(outdir + '/prob_window_variability_contig' + str(contig) + '.txt', 'a+')
		f_mc = open(outdir + '/prob_main_chain_contig' + str(contig) + '.txt', 'a+')

		for gene in data[0]:
			f_wc.write(self.genes_decode[gene] + '\t' + str(data[0][gene]) + '\n')
			f_mc.write(self.genes_decode[gene] + '\n')

		for gene in data[1]:
			f_io.write(self.genes_decode[gene] + '\t' + str(data[1][gene]) + '\n')

		for bridge in data[2]:
			f_ab.write(self.genes_decode[bridge[0]] + '\t' + self.genes_decode[bridge[1]] + '\t' + str(data[2][bridge]) + '\t')


	def compute_variability(self, outdir, reference, window=20, iterations=500, min_depth=0, max_depth=-1):
		
		print('Reference is ' + reference)
		print('Number of contigs: ' + str(len(self.list_graph[reference])))

		for contig in range(len(self.list_graph[reference])):
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
						for i in range(min(start_index, end_index) + 1, max(start_index, end_index)):
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
						for i in range(min(start_index, end_index) + 1, max(start_index, end_index)):
							prob_variability_table[base_line[i]] += 1/float(norm)

				print('\r', end='')
				count += 1

			self.save_data([variability_table, io_table, all_bridges, 
						prob_variability_table, prob_io_table, prob_all_bridges], outdir, contig)
	
		print('\nComputing completed')