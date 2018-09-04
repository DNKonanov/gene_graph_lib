import GC_skew
from numpy import trapz
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--base_complexity', type=str, default='no', help='file with base ori aligned complexity')
parser.add_argument('--fna_file', type=str, default='no', help='full fna file')
parser.add_argument('--out_file', type=str, default='no', help='output file with GC aligned complexity.')
parser.add_argument('--reverse', type=int, default=0, help='if 1, out complexity table will be reversed. Default 0')

args = parser.parse_args()

try:
	base_complexity_file = open(args.base_complexity, 'r')
	fna_file = open(args.fna_file, 'r')

except FileNotFoundError:
	print('Incorrect input!')



orient = int(args.reverse)
base = GC_skew.get_base(args.fna_file)
tmp_x = []
for line in base_complexity_file:
	tmp_x.append(float(line.split('\t')[1][:-1]))

k = len(tmp_x)/10000.0

x = []
count = 0
while(len(x) < 10000 and count <= len(tmp_x)):
	x.append(max(tmp_x[int(count):int(count+k)]))
	count += k



if orient == 1:
	x.reverse()
	base = base[::-1]

x = x/trapz(x, x=base)


f_out = open(args.out_file, 'a+')
for i in range(len(x)):
	f_out.write(str(base[i]) + '\t' + str(x[i]) + '\n')


