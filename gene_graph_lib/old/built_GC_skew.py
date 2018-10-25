import GC_skew
import matplotlib.pyplot as plt
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument('--fna', type=str, default='no', help='full fna file')

args = parser.parse_args()

for rec in SeqIO.parse(args.fna, 'fasta'):
        
    gc = GC_skew.skew_count(rec)

    ori = gc.index(min(gc))
    break

print('Origin coordinate: ' + str(ori))

plt.xlabel('bp')
plt.ylabel('GC skew')
plt.plot(gc)
plt.title('GC skew')
plt.axvline(x=ori, c='r', label='Origin')
plt.legend()
plt.show()

