import random
import numpy as np
from numpy import trapz
from Bio import SeqIO


def get_base(fna_file):
    for rec in SeqIO.parse(fna_file, 'fasta'):
        
        gc = skew_count(rec)
        
        
        last = gc[-1]
        N = len(gc)
        k = N/10000
        new_gc = []
        i = 0
        while i <= N - 1:
            new_gc.append(gc[int(i)] - last*(i/N))
            i += k


        ori = new_gc.index(min(new_gc))
        
        new_gc = new_gc[ori:] + new_gc[:ori]
        
        
        r = new_gc.index(max(new_gc))/len(new_gc)
        skew = 0.5 - r
        
        new_gc = np.array(new_gc) - min(new_gc)
        new_gc = new_gc/max(new_gc)     
        
        base = np.array([i/len(new_gc) for i in range(len(new_gc))])
        
        for coord in range(len(base)):
            if base[coord] <= r:
                base[coord] += skew*base[coord]/r
            
            elif base[coord] > r:
                base[coord] += skew*((1 - base[coord])/(1 - r))
        
        return base
        break
                
def skew_count(inseq):
    Skew_counter = []    
    skew = 0
    Skew_counter.append(skew)
    for i in inseq:
        if i == 'C' or i == 'c':
            skew -= 1
            Skew_counter.append(skew)
        elif i == 'G' or i == 'g':
            skew += 1
            Skew_counter.append(skew)
        else:
            Skew_counter.append(skew)
    return Skew_counter