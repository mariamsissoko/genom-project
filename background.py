
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 13:58:47 2016

@authors: Maria Virginia Ruiz & Mariam Sissoko 
"""
import matplotlib.pyplot as plt
import numpy as np
import os
path_dir = os.getcwd() 
import random
import dicNames
from igraph import *
import scipy.stats as scp
import scipy.spatial as scp2 
import scipy
import math


'''
Permet de lire un fichier contenant un genome au format fasta 
#Input :nom du fichier fasta 
#Output:chaine de caractere contenant la sequence 
'''
mots = [ 'N','Y','R', 'W', 'S', 'B', 'V', 'K', 'M', 'H', 'D']

def lireSeq(infile): 
    f = open(infile, 'r')
    name = ""
    seq = ""
    for line in f:
        if ">" not in line:
            seq =seq+line[0:len(line)-1]
        else:
            if "Mycoplasma genitalium" not in line and "Pyrococcus abyssi GE5" not in line:
                name = line[12:len(line)-1]
            else:
                if "Mycoplasma genitalium" not in line:
                    name = line[10:len(line)-1]
                else:
                    name = line[16:len(line)-1]
    f.close()
    for i in mots:
        if i in seq:
            seq = seq.replace(i, '')
    nameToReturn = dicNames.dicNames[(name.split(',')[0])]+': '+name.split(',')[0]
    return nameToReturn, seq

#name, seq = lireSeq(path_dir+"\\GCA_000008085.1_ASM808v1_genomic.fna")

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
    kmers = dict(dico)
    n = len(string)
    for i in range(0, n-k+1):
        key = string[i:i+k]
        kmers[key] += 1
    return kmers


    
'''
Proportion de kmers dans un genome donne
'''
def proportion_of_kmers(nb_kmers, dic_kmers):
    kmers_prop = dict(dic_kmers)
    for i in dic_kmers.keys():
        kmers_prop[i]/=float(nb_kmers)
    return kmers_prop

def proportion_along_genome(seq, pas, fenetre, dico_init_kmers, k):
    liste_prop=[]
    for i in range(0,len(seq)-fenetre,pas):
        dico = dict(dico_init_kmers) 
        dico = find_kmers(dico, k,seq[i:i+fenetre] )
        dico = proportion_of_kmers(len(seq[i:i+fenetre])-k+1,dico)
        liste_prop.append(dico)
    return liste_prop

def distanceEuclidienne(genome, partieGenome):
    dist = 0.0
    for i in genome.keys():
        dist = dist + ((genome[i]-partieGenome[i])**2)
    return dist


def do_mean_sub_parts(dic_sub_parts):
    dic_mean={}
    for i in dic_sub_parts.keys():
        dic_mean[i]=sum(dic_sub_parts[i])/len(dic_sub_parts[i])
    return dic_mean 

def do_sd_sub_parts(dic_sub_parts, dic_means):
    dic_sd={}
    for i in dic_sub_parts.keys():
        dic_sd[i] = 0
        for j in range(len(dic_sub_parts[i])):
            dic_sd[i] = dic_sd[i] + (dic_sub_parts[i][j]-dic_means[i])**2
            
        dic_sd[i]= math.sqrt(dic_sd[i]/(len(dic_sub_parts[i])-1))
    return dic_sd

def calculeDistanceParties(genome, partiesGenome):
    distances = []
    for i in range(len(partiesGenome)):
        distances.append(distanceEuclidienne(genome, partiesGenome[i]))
    return distances

def calculeDistancePartiesCorr(genome, partiesGenome):
    distances = []
    for i in range(len(partiesGenome)):
        x = scipy.array(list(genome.values()))
        y = scipy.array(list(partiesGenome[i].values()))
        distances.append(scp2.distance.correlation(x, y))
    return distances

def calculeDistancePartiesCorr2(genome, partieGenome):
    distance = 0.0
    x = scipy.array(list(genome.values()))
    y = scipy.array(list(partieGenome.values()))
    distance = scp2.distance.correlation(x, y)
    return distance
        
def do_distances_all_genomes(genomes, partiesGenomes):
    distances = {}
    for i in genomes.keys():
        distances[i] = distanceEuclidienne(genomes[i], partiesGenomes[i])
    return distances
    
def do_profils_all_genomes(dic_genomes, dico_k, k):
    dic_profils={}
    for i in dic_genomes.keys():
        dic_kmers = find_kmers(dico_k, k, dic_genomes[i])
        dic_profils[i] = proportion_of_kmers(len(dic_genomes[i])-k+1, dic_kmers)
    return dic_profils 

def do_profil_genome(genome, dico_k, k):
    dic_kmers = find_kmers(dico_k, k, genome)
    dic_profil = proportion_of_kmers(len(genome)-k+1, dic_kmers)
    return dic_profil 

def do_sub_parts(profils, dic_genomes, pas, fenetre, dico_init_kmers, k , euclidienne):
    dic_sub_parts={}
    for i in profils.keys():
        prop = proportion_along_genome(dic_genomes[i], pas, fenetre, dico_init_kmers, k)
        if euclidienne == 1:
            dic_sub_parts[i] = calculeDistanceParties(profils[i], prop)
        else:
            dic_sub_parts[i] = calculeDistancePartiesCorr(profils[i], prop)
        
    return dic_sub_parts

def get_partieDif_all_genomes(genomes, positions, pas, fenetre, cont):
    partieDif = {}
    try:
        for i in genomes.keys():
            seq = genomes[i]
            position = positions[i]
            posIni = int(position[cont] - (fenetre/2))
            if  posIni < 0 :
                posIni = 0
            posFini = int(position[cont] + (fenetre/2))
            if posFini > len(seq):
                posFini = len(seq)
            else: 
                if posFini < 0 :
                    posFini = fenetre
            #print (i, posIni, posFini)
            partieDif[i] = seq[posIni : posFini ]
    except Exception as e:
        print ("ici 5", e)
        
    return partieDif


def transfert_study(dic_profils, dic_parts_genomes, dic_kmers, k, euclidienne):
    dic_profils_part = do_profils_all_genomes(dic_parts_genomes, dic_kmers, k)
    matrice = np.eye(len (dic_profils_part.keys()), len (dic_profils_part.keys()))
    keys = dic_profils_part.keys()
    k=0
    l=0
    for i in keys:
        for j in keys:
            if i!=j:
                if euclidienne == 1:
                    matrice[k,l] = distanceEuclidienne(dic_profils[j], dic_profils_part[i])
                else:
                    matrice[k,l] = calculeDistancePartiesCorr2(dic_profils[j], dic_profils_part[i])
            l+=1
        k+=1
        l=0
    return matrice, list(keys) 

def tranfert_study_result(matrix, cles, infoDistances):
    liste_probable_transfert = []
    for i in range(len(matrix)):
        transfert = True
        genome = cles[i]
        
        info = infoDistances[genome] 
        for j in info:
            if j[0] < j[1]:
                transfert = False
                break
        if transfert:
            genome2 = cles[np.argmin(matrix[i])]
            liste_probable_transfert.append([genome, genome2])
    return liste_probable_transfert

def transfert_studyBoutGenome(dic_profils, butGenome, euclidienne):
    liste = []
    keys = list(dic_profils.keys())
    for i in dic_profils.keys():
        if euclidienne == 1:
            liste.append(calculeDistanceParties(dic_profils[i], butGenome))
        else:
            liste.append(calculeDistancePartiesCorr2(dic_profils[i], butGenome))
    indx = np.argmin(liste)
    genomePredicte =  keys[indx]
    return genomePredicte

def predict_Genome(dic_profils, sequences, dic_kmers, k, taille, euclidienne, typePredict=False):
    listDist=[]
    if not typePredict:
        genome = random.choice(list(dic_profils.keys()))
        posGenome = random.randint(0,len(sequences[genome]) - taille)
        boutGenome = [genome, sequences[genome][posGenome:posGenome+taille]]
    else:
        posGenome = random.randint(0,len(sequences) - taille)
        boutGenome = ['', sequences[posGenome:posGenome+taille]]
    
    profilMorceau = do_profil_genome(boutGenome[1], dic_kmers, k)
            
    keys = list(dic_profils.keys())
    for j in keys:
        if euclidienne == 1:
            listDist.append(distanceEuclidienne(dic_profils[j], profilMorceau))
        else :
            listDist.append(calculeDistancePartiesCorr2(dic_profils[j], profilMorceau))
    boutGenome.append(keys[listDist.index(min(listDist))])
    if not typePredict:
        return boutGenome, genome
    else:
        return boutGenome

def experiment(dic_profils, sequences, dic_kmers, k, exp, taille, euclidienne):
    vectResult=[]
    for i in range(exp):
        res, genome = predict_Genome(dic_profils, sequences, dic_kmers, k, taille, euclidienne)
        if res[0] == res[2]:
            vectResult.append(1)
        else:
            vectResult.append(0)
    return vectResult, genome

def errorRate(vectResult):
    rate=0
    for i in range(len(vectResult)):
        if vectResult[i]==0:
            rate+=1
    return rate/float(len(vectResult))

def experiment2(dic_profils, sequences, dic_kmers, k, exp, taille, euclidienne):
    vectResult = []
    tauxDerreur = []
    keys = list(dic_profils.keys())
    for j in keys:
        resulst = []
        for i in range(exp):
            res = predict_Genome(dic_profils, sequences[j], dic_kmers, k, taille, euclidienne, True)
            if j == res[2]:
                resulst.append(1)
            else:
                resulst.append(0)
        vectResult.append(resulst)
    for i in vectResult:
        tauxDerreur.append(errorRate(i))
    return tauxDerreur, keys

def sequencesDif(liste):
    f=open("SeqDifs.txt",'w')
    for i in liste.keys():
        f.write(i+'\n')
        f.write(liste[i]+'\n')
    f.close()

def file_result(liste_transfert, file_name):
    f=open(file_name,'w')
    relations = []
    f.write("Result for transfert analysis\n")
    for i in liste_transfert:
        f.write('For [ '+i[0]+' ] probable transfert from [ '+i[1]+' ]\n')
        relations.append([i[0], i[1]])
    f.close()
    
    return relations
    
def do_data_matrix(profils_genomes):
    liste = list(profils_genomes.keys())
    data_matrix = np.zeros((len(profils_genomes), len(profils_genomes[liste[0]])))
    for i in liste:
        profil = profils_genomes[i]
        liste2 = list(profil.keys())
        for j in liste2:
            data_matrix[liste.index(i),liste2.index(j)] = profil[j]
    listToReturn = []
    for i in liste:
        typeS = i.split(':')[0]
        listToReturn.append(typeS)
    return data_matrix, listToReturn

def reseau2(relations, listeNodes, types):
    g = Graph(directed=True)
    g.add_vertices(listeNodes)
    #g.vs['label'] = types
    for i in relations:
        g.add_edges([(i[0],i[1])])
    
    layout = g.layout("fruchterman_reingold")
    color_dict = {"Archaea": "blue", "Bacteria": "pink"}
    g.vs["color"] = [color_dict[t] for t in types]
    plot(g, "Reseau.pdf", layout= layout)