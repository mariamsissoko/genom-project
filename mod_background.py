



# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 13:58:47 2016

@authors: Maria Virginia Ruiz & Mariam Sissoko 
"""
import matplotlib.pyplot as plt
import os
import numpy as np 
import rpy 
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
def find_kmers(dico, k, string ):
    kmers =dict(dico)
    n = len(string)
    for i in range(0, n-k+1):
        key = string[i:i+k]
        kmers[key] += 1
    return kmers

'''
Proportion de kmers dans un genome donne
Input :nombre de kmers dans le genome,dictionnaire de comptage ds kmers 
Output:Dictionnare contenant la proportion de chaque kmer 
'''
def proportion_of_kmers(nb_kmers, dic_kmers):
    kmers_prop=dict(dic_kmers)
    for i in dic_kmers.keys():
        kmers_prop[i]/=float(nb_kmers)
    return kmers_prop

'''
Fenetre glissante pour calculer la propotion de kmers dans chaque partie du genome

'''
def proportion_along_genome(seq,pas,fenetre,dico_init_kmers,k):
    liste_prop=[]
    for i in range(0,len(seq)-fenetre,pas):
        dico=dict(dico_init_kmers) 
        dico=find_kmers(dico, k,seq[i:i+fenetre] )
        dico=proportion_of_kmers(len(seq[i:i+fenetre])-k+1,dico)
        liste_prop.append(dico)
    return liste_prop


def distanceEuclidienne(genome, partieGenome):
    dist = 0.0
    for i in genome.keys():
        dist=dist+((genome[i]-partieGenome[i])**2)
    return dist


def calculeDistanceParties(genome, partiesGenome):
    distances = []
    for i in range(len(partiesGenome)):
        distances.append(distanceEuclidienne(genome, partiesGenome[i]))
    return distances

def do_profils_all_genomes(dic_genomes,dico_k,k):
    dic_profils={}
    for i in dic_genomes.keys():
        dic_kmers=find_kmers(dico_k, k,dic_genomes[i])
        dic_profils[i]=proportion_of_kmers(len(dic_genomes[i])-k+1,dic_kmers)
    return dic_profils 
        
def do_distance_matrix(profils):
    mat=np.array(len(profils),len(profils))
    k=0
    l=0
    for i in profils.keys():
        for j in profils.keys():
            mat[k][l]=distanceEuclidienne(profils[i],profils[j])
            l+=1
        k+=1
    return mat 

def do_sub_parts(profils,dic_genomes,pas,fenetre,dico_init_kmers,k):
    dic_sub_parts={}
    for i in profils.keys():
        prop=proportion_along_genome(seq,pas,fenetre,dico_init_kmers,k)
        dic_sub_parts[i]=calculeDistanceParties(genome,prop)
    return dic_sub_parts

def plot_for_each_genome(sub_parts_for_each_genome):
    fig = plt.figure()
    number1=np.ceil(np.sqrt(len(sub_parts_for_each_genome)))
    number2=len(sub_parts_for_each_genome)
    plt.subplot(number1number1number2)
    for i in sub_parts_for_each_genome.keys():
        x_points = range(0,len(sub_parts_for_each_genome[i]))
        y_points =sub_parts_for_each_genome[i]
        p = ax.plot(x_points, y_points, 'b')
        ax.set_xlabel('position ')
        ax.set_ylabel('distance')
        #ax.set_title('distances entre le profil moy')
        fig.show()
        
def transfert_study(dic_profils,dic_parts_genomes,dic_kmers,k):
    dic_profils_part = do_profils_all_genomes(dic_genomes,dic_kmers,k)
    matrice = np.eyes(len (dic_profils_part.keys()), len (dic_profils_part.keys()))
    listCles = {}
    cont = 0
    count+=1
    keys = dic_profils_part.keys()
    for i in keys:
        listCles[i] = cont
    for i in keys:
        for j in keys:
            if i!=j:
                matrice = distanceEuclidienne(genome[j], partieGenome[i])
    return matrice 

def tranfert_study_result(list_genomes_names,matrix):
    liste_probable_transfert=[]
    for i in range(len(matrix)):
        liste_probable_transfert.append(np.argmin(matrix(i)))
    return liste_probable_transfert 

def read_file(file_name):
    f=open(file_name)
    ligne=f.readline()
    dic_type={}
    while len(ligne)>0:
        ligne=ligne.split(' ')
        name=""
        for i in range(len(ligne)-1):
            name=name+ligne[i]+" "
        name=name[0:len(name)-1]
        dic_type[name]=ligne[len(ligne)-1][0:len(ligne[len(ligne)-1])-1]
        ligne=f.readline()
    return dic_type

def file_result(liste_transfert,liste_genomes_names,dic_type,file_name):
    f=open(file_name,'w')
    f.write("result for transfert analysis\n")
    for i in range(len(liste_transfert)):
        f.write('for '+liste_genomes_names[i]+' probable tranfert from '+liste_genomes_names[liste_transfert[i]]+'\n')
        f.write('transfert from '+dic_type[liste_genomes_names[i]]+' to '+dic_type[liste_genomes_names[liste_transfert[i]]]+'\n\n')
    return 

def do_data_matrix(profils_genomes,nb_kmers)
    data_matrix=np.zeros(len(profils_genomes),nb_kmers+1)
    liste=profil_genomes.keys()
    for i in profils_genomes.keys():
        data_matrix[liste.index(i),0]=i
        for j in range(1,len(profils_genomes[i])+1):
             data_matrix[liste.index(i),j]=profils_genomes[i][j]
    return data_matrix 
            
    
data_matrix=do_data_matrix()
distance_matrix=r.dist(as.matrix(data_matrix))
genome_hclust=r.hclust(distance_matrix)
r.plot(genom_hclust)

"""
#print (path_dir)
k=6
seq = lireSeq(path_dir+"\\seq.fasta")
alphabet_nt={'A':0,'T':0,'G':0,'C':0}

dico_k_mers=do_dic_kmers(k,alphabet_nt)
kmers = find_kmers(dico_k_mers,k,seq)
#print (kmers)
genome = proportion_of_kmers(len(seq)-k+1, kmers)

parties = proportion_along_genome(seq,10000,100000,dico_k_mers,k)

#print (genome)
print (len(parties))

distances = calculeDistanceParties(genome, parties)

fig = plt.figure()
ax = fig.add_subplot(111)
x_points = range(0,len(parties))
y_points = distances
p = ax.plot(x_points, y_points, 'b')
ax.set_xlabel('x-points')
ax.set_ylabel('y-points')
ax.set_title('Simple XY point plot')
fig.show()

position = distances.index(max(distances))
print (position)

partieGenomeDif = seq[position*10000: position*10000+100000]
#print (len(partieGenomeDif))

#f = open(infile, 'w')
subParties = proportion_along_genome(partieGenomeDif,100,1000,dico_k_mers,k)
distances2 = calculeDistanceParties(genome,subParties)

fig = plt.figure()
ax = fig.add_subplot(111)
x_points = range(0,len(subParties))
y_points = distances2
p = ax.plot(x_points, y_points, 'b')
ax.set_xlabel('x-points')
ax.set_ylabel('y-points')
ax.set_title('Simple XY point plot')
fig.show()

position2 = distances2.index(max(distances2))
print (position2)

pos = (position*10000) + (position2*100)

print (pos)

partieGenomeDif = seq[pos: pos+1000]

len(partieGenomeDif)

subParties = proportion_along_genome(partieGenomeDif,50,100,dico_k_mers,k)
distances2 = calculeDistanceParties(genome,subParties)

fig = plt.figure()
ax = fig.add_subplot(111)
x_points = range(0,len(subParties))
y_points = distances2
p = ax.plot(x_points, y_points, 'b')
ax.set_xlabel('x-points')
ax.set_ylabel('y-points')
ax.set_title('Simple XY point plot')
fig.show()

position3 = distances2.index(max(distances2))
print (position3)

pos = (position*10000) + (position2*100) + (position3*50)

print (pos)

partieGenomeDif = seq[pos: pos+100]
len(partieGenomeDif)

"""
