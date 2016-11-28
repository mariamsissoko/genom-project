# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 13:58:47 2016

@authors: Maria Virginia Ruiz & Mariam Sissoko 
"""
import os
path_dir = os.getcwd()  

'''
Permet de lire un fichier contenant un genome au format fasta 
#Input :nom du fichier fasta 
#Output:chaine de caractere contenant la sequence 
'''
def lireSeq(infile):
    f = open(infile, 'r')
    seq = ""
    for line in f:
        if ">" not in line:
            seq =seq+line[0:len(line)-1]
    f.close()
    return seq

'''
Permet de faire la liste de tous les kmers possibles pour un k donne
et initialise le comptage à 0 
Input :longueur des kmers ,alphabet (longueur 4 pour de l'ADN)
Output:dictionnaire avec les kmers possibles comme clés et 0 comme valeur
'''
def do_dic_kmers(k,alphabet):
    dic_k_mers=alphabet
    for k in range(k-1):
        new_dic_k_mers={}
        for i in dic_k_mers.keys():
            for j in alphabet.keys():
                new_dic_k_mers[i+j]=0
        dic_k_mers=new_dic_k_mers
    return dic_k_mers

'''
Comptage des kmers d'une sequence donnée
Input :Dictionnaire avec tous les kmers possibles, k , sequence d'ADN
Output:Dictionnaire qui a pour clé un kmer et comme clé associé le nombre de fois ou le kmer apparait
'''
def find_kmers(dico, k,string ):
    kmers =dict(dico)
    n = len(string)
    for i in range(0, n-k+1):
        key = string[i:i+k]
        kmers[key] += 1
    return kmers

'''
Proportion de kmers dans un genome donne
'''
def proportion_of_kmers(nb_kmers,dic_kmers):
    kmers_prop=dict(dic_kmers)
    for i in dic_kmers.keys():
        kmers_prop[i]/=float(nb_kmers)
    return kmers_prop

def proportion_along_genome(seq,pas,fenetre,dico_init_kmers,k):
    liste_prop=[]
    for i in range(0,len(seq)-fenetre,pas):
        dico=dict(dico_init_kmers) 
        dico=find_kmers(dico, k,seq[i:i+fenetre] )
        dico=proportion_of_kmers(len(seq[i:i+fenetre])-k+1,dico)
        liste_prop.append(dico)
    return liste_prop
        
print (path_dir)
k=2
seq = lireSeq(path_dir+"\\seq.fasta")
alphabet_nt={'A':0,'T':0,'G':0,'C':0}
dico_k_mers=do_dic_kmers(k,alphabet_nt)
kmers = find_kmers(dico_k_mers,k,seq)
print kmers
kmers_prop=proportion_of_kmers(len(seq)-k+1,kmers)
print kmers_prop
print len(proportion_along_genome(seq,10000,100000,dico_k_mers,k))
