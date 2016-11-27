# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 13:58:47 2016

@author: MaVi
"""
import os
path_dir = os.getcwd()  

def lireSeq(infile):
    f = open(infile, 'r')
    seq = ""
    for line in f:
        if ">" not in line:
            seq = line
    f.close()
    return seq
    
print (path_dir)
seq = lireSeq(path_dir+"\\seq.fasta")

def find_kmers(string, k):
    kmers = {}
    n = len(string)
    for i in range(0, n-k+1):
        key = string[i:i+k]
        if key in kmers:
            kmers[key] += 1
        else: 
            kmers[key] = 1
    return kmers

kmers = find_kmers(seq,2)

