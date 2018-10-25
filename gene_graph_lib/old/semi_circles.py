#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 12:59:31 2017

@author: dmitry
"""

import cairo
import math
import argparse

parser = argparse.ArgumentParser(prog='create_picture.sh')
parser.add_argument('-i', '--input_file', default='no', type=str, help='input_file')
parser.add_argument('-o', '--output_dir', default='no', type=str, help='ounput directory')
parser.add_argument('--freq_min', type=int, default=0, help='Minimal count of paths for drawing')
parser.add_argument('--complexity_file', type=str, default='no', help='Complexity file (additive or window)')
args = parser.parse_args()


def draw_semicircles_bonds(freq_table, main_chain, freq_min=0, complexity=[], selected_genes=[]):

	surface = cairo.SVGSurface(args.output_dir+'/graph_' + args.input_file.split('/')[-1].split('.')[0] + '.svg', (len(main_chain) + 1)*10, 10*len(main_chain))
	context = cairo.Context(surface)
	X_coord = []
	Y_coord = []
	Xi_Xj_bonds = []
	dict_freq_table = {}

	print(freq_table)

	for i in freq_table:
		dict_freq_table.update([(i[0]+i[1], i[2])])
	print(dict_freq_table)
	max_freq = 0
	for i in freq_table:
		if i[2]  > max_freq:
			max_freq = i[2]
	print(max_freq)
	x = 10
	y = 10*len(main_chain)/2
	for gene in main_chain:
		print(gene)
		X_coord.append(x)
		x += 10
		Y_coord.append(y)
	N = len(freq_table)
	count = 1
	#draw semicircles
	context.set_line_width(1.6)
	for i in range(len(main_chain)):
		for j in range(i + 1 , len(main_chain)):
			if main_chain[i]+main_chain[j] in dict_freq_table:
				if dict_freq_table[main_chain[i]+main_chain[j]] >= freq_min:
					
					print(str(count) + ' from ' + str(N))
					count += 1
					context.set_source_rgba(0, 0, 1, ((0.75)*dict_freq_table[main_chain[i]+main_chain[j]]/float(max_freq) + 0.25))
					context.arc_negative((X_coord[j] + X_coord[i])/2, Y_coord[0], (X_coord[j] - X_coord[i])/2, 0, math.pi)
					print('Added')
					context.stroke()
			elif main_chain[j]+main_chain[i] in dict_freq_table:
				if dict_freq_table[main_chain[j]+main_chain[i]] >= freq_min:
					
					print(str(count) + ' from ' + str(N))
					count += 1
					context.set_source_rgba(0, 0, 1, ((0.75)*dict_freq_table[main_chain[j]+main_chain[i]]/float(max_freq) + 0.25))
					context.arc((X_coord[j] + X_coord[i])/2, Y_coord[0], (X_coord[j] - X_coord[i])/2, 0, math.pi)
					print('Added')
					context.stroke()


	context.move_to(X_coord[0], Y_coord[0])

	#draw complexity
	if len(complexity) > 0:
		context.set_line_width(10)
		context.set_source_rgba(0,1,0,1)
		if complexity != []:
			for i in range(len(complexity)):
				print(X_coord[i], Y_coord[0] - 10*complexity[i][1])
				context.line_to(X_coord[i], Y_coord[0] - 10*complexity[i][1])
			context.stroke()


	#draw main chain
	context.set_line_width(1.6)
	context.set_source_rgba(0, 0, 0, 1)
	for i in range(len(X_coord)):
		context.arc(X_coord[i], Y_coord[i], 0.8, 0, 2 * math.pi)
		context.stroke()
	context.move_to(X_coord[0], Y_coord[0])
	context.line_to(X_coord[-1], Y_coord[-1])
	context.stroke()


	
	#append stamms
	FOUND = 0
	context.set_font_size(1)
	for i in range(len(main_chain)):
		if main_chain[i] in selected_genes:
			FOUND += 1
			context.set_line_width(10)
			context.move_to(X_coord[i], Y_coord[i] + 1000)
			context.set_source_rgba(1,0,0,1)
			context.line_to(X_coord[i], Y_coord[i] + 50)
			context.stroke()
			context.set_source_rgba(0,0,0,1)
		context.move_to(X_coord[i] - 3, Y_coord[i] + 3)
		context.show_text(main_chain[i])

	
	print FOUND

	

	surface.finish()

def main(input_freq_table, input_complexity_file, freq_min):
	if input_freq_table == 'no':
		print 'NO INPUT FILE'
		return 0

	f_in = open(input_freq_table)
	freq_table = [[], []] 
	for line in f_in:

		string = line.split('\t')
		print(string)
		if len(string) == 3:
			
	
			gene1 = string[0]
	
			gene2 = string[1]
			freq = int(string[2])

			freq_table[1].append([gene1, gene2, freq])
		

		elif line != '':
			print(line)
			freq_table[0].append(line[:-1])
	selected_genes = [] #will be added
	complexity = []
	if input_complexity_file != 'no':
		file_complexity = open(input_complexity_file, 'r')

		
		for line in file_complexity:
			if line.count('\t'):
				complexity.append([line.split('\t')[0], int(round(float(line.split('\t')[1][:-1])))])
	draw_semicircles_bonds(freq_table[1], freq_table[0], freq_min=freq_min, complexity=complexity)



if args.input_file != 'no':
	main(args.input_file, args.complexity_file, args.freq_min)