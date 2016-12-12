
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
#import pygraphviz as pgv
from igraph import *
import scipy.stats as scp


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
    print (name.split(',')[0])
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


def calculeDistanceParties(genome, partiesGenome):
    distances = []
    for i in range(len(partiesGenome)):
        distances.append(distanceEuclidienne(genome, partiesGenome[i]))
    return distances

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

def do_sub_parts(profils, dic_genomes, pas, fenetre, dico_init_kmers, k):
    dic_sub_parts={}
    for i in profils.keys():
        prop = proportion_along_genome(dic_genomes[i], pas, fenetre, dico_init_kmers, k)
        dic_sub_parts[i] = calculeDistanceParties(profils[i], prop)
        print (scp.mannwhitneyu(profils[i], prop, use_continuity=True))
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


def transfert_study(dic_profils, dic_parts_genomes, dic_kmers, k):
    dic_profils_part = do_profils_all_genomes(dic_parts_genomes, dic_kmers, k)
    matrice = np.eye(len (dic_profils_part.keys()), len (dic_profils_part.keys()))
    keys = dic_profils_part.keys()
    k=0
    l=0
    for i in keys:
        for j in keys:
            if i!=j:
                matrice[k,l] = distanceEuclidienne(dic_profils[j], dic_profils_part[i])
            l+=1
        k+=1
        l=0
    return matrice, list(keys) 

def tranfert_study_result(matrix):
    liste_probable_transfert = []
    for i in range(len(matrix)):
        liste_probable_transfert.append(np.argmin(matrix[i]))
    return liste_probable_transfert 

def transfert_studyBoutGenome(dic_profils, butGenome):
    liste = []
    keys = list(dic_profils.keys())
    for i in keys:
        liste.append(distanceEuclidienne(dic_profils[i], butGenome))
    indx = np.argmin(liste)
    genomePredicte =  keys[indx]
    return genomePredicte

def predict_Genome(dic_profils, sequences, dic_kmers, k, taille, typePredict=False):
    listDist=[]
    if not typePredict:
        genome = random.choice(list(dic_profils.keys()))
        posGenome = random.randint(0,len(sequences[genome]) - taille)
        boutGenome = [genome, sequences[genome][posGenome:posGenome+taille]]
    else:
        posGenome = random.randint(0,len(sequences) - taille)
        boutGenome = ['', sequences[posGenome:posGenome+taille]]
    
    profil_boutGenome = find_kmers(dic_kmers, k, boutGenome)
    profil_boutGenome = proportion_of_kmers(len(boutGenome[1])-k+1, dic_kmers)
    keys = list(dic_profils.keys())
    for j in keys:
        listDist.append(distanceEuclidienne(dic_profils[j],profil_boutGenome))
    boutGenome.append(keys[listDist.index(min(listDist))])
    if not typePredict:
        return boutGenome, genome
    else:
        return boutGenome

def experiment(dic_profils, sequences, dic_kmers, k, exp, taille):
    vectResult=[]
    for i in range(exp):
        res, genome = predict_Genome(dic_profils, sequences, dic_kmers, k, taille)
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

def experiment2(dic_profils, sequences, dic_kmers, k, exp, taille):
    vectResult = []
    tauxDerreur = []
    keys = list(dic_profils.keys())
    for j in keys:
        resulst = []
        for i in range(exp):
            res = predict_Genome(dic_profils, sequences[j], dic_kmers, k, taille, True)
            if j != res[2]:
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

def file_result(liste_transfert, liste_genomes_names, file_name):
    f=open(file_name,'w')
    relations = []
    f.write("Result for transfert analysis\n")
    for i in range(len(liste_transfert)):
        f.write('For [ '+liste_genomes_names[i]+' ] probable transfert from [ '+liste_genomes_names[liste_transfert[i]]+' ]\n')
        relations.append([liste_genomes_names[i], liste_genomes_names[liste_transfert[i]]])
    f.close()
    
    return relations
    
def do_distance_matrix(profils):
    matrice = np.zeros((len(profils.keys()),len(profils.keys())))
    k=0
    l=0
    keys = list(profils.keys())
    for i in keys:
        for j in keys:
            matrice[k, l] = distanceEuclidienne(profils[i], profils[j])
            l+=1
        k+=1
        l=0
    
    return matrice

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


def reseau(relations, listeNodes):
    graph = pgv.AGraph(directed=True)
    # changement de l'attribut shape (forme) pour le type box (boîte)
    #graph.node_attr["shape"] = "box"
    # le premier noeud doit être rempli (style="filled") et coloré en gris (color='grey')
    """
    for i in listeNodes:
        
        if i == "Archaea":
            graph.add_node(i, style="filled", color='green')
        else: 
            graph.add_node(i, style="filled", color='blue')
    """
    for i in relations:
        graph.add_edge(i[0], i[1])
    graph.layout('dot')
    graph.draw("/images/Reseaux.pdf")
    graph.close()
    

def reseau2(relations, listeNodes, types):
    g = Graph(directed=True)
    g.add_vertices(listeNodes)
    #g.vs['label'] = types
    
    for i in relations:
        g.add_edges([(i[0],i[1])])
    
    layout = g.layout_kamada_kawai()
    color_dict = {"Archaea": "blue", "Bacteria": "pink"}
    g.vs["color"] = [color_dict[t] for t in types]
    plot(g, "Reseau.pdf", layout= layout)

relations = [['Archaea: Desulfurococcus kamchatkensis 1221n', 'Bacteria: Thermobaculum terrenum ATCC BAA-798 chromosome 2'], ['Bacteria: Shewanella putrefaciens CN-32', 'Bacteria: Bacillus pseudofirmus OF4 plasmid pBpOF4-02'], ['Archaea: Methanosarcina mazei strain Goe1', 'Bacteria: Acholeplasma laidlawii PG-8A'], ['Archaea: Natronomonas pharaonis DSM 2160 plasmid PL23 complete genome', 'Bacteria: Chloroflexus aurantiacus J-10-fl'], ['Archaea: Methanocorpusculum labreanum Z', 'Bacteria: Bdellovibrio bacteriovorus complete genome'], ['Bacteria: Methylacidiphilum infernorum V4', 'Bacteria: Carboxydothermus hydrogenoformans Z-2901'], ['Bacteria: Acidobacterium capsulatum ATCC 51196', 'Bacteria: Opitutus terrae PB90-1'], ['Archaea: Pyrobaculum aerophilum str. IM2', 'Archaea: Pyrobaculum arsenaticum DSM 13514'], ['Archaea: Methanopyrus kandleri AV19', 'Archaea: Caldivirga maquilingensis IC-167'], ['Archaea: Methanoculleus marisnigri JR1', 'Archaea: Methanocorpusculum labreanum Z'], ['Bacteria: Chlamydia trachomatis strain L2/434/Bu complete genome', 'Archaea: Pyrococcus furiosus DSM 3638'], ['Bacteria: Rickettsia rickettsii str. Iowa', 'Bacteria: Acholeplasma laidlawii PG-8A'], ['Bacteria: Chlorobium chlorochromatii CaD3', 'Bacteria: Shewanella putrefaciens CN-32'], ['Bacteria: Nostoc punctiforme PCC 73102 plasmid pNPUN05', 'Bacteria: Shewanella putrefaciens CN-32'], ['Archaea: Haloquadratum walsbyi DSM 16790 plasmid PL47 complete genome', 'Bacteria: Corynebacterium glutamicum ATCC 13032 DNA'], ['Archaea: Metallosphaera sedula DSM 5348', 'Archaea: Archaeoglobus profundus DSM 5631 plasmid pArcpr01'], ['Bacteria: Dictyoglomus thermophilum H-6-12', 'Bacteria: Carboxydothermus hydrogenoformans Z-2901'], ['Archaea: Methanococcus maripaludis C6', 'Bacteria: Bacillus cereus Q1 plasmid pBc53'], ['Bacteria: Mycoplasma genitalium G37 strain G-37 unique_12', 'Bacteria: Borrelia burgdorferi B31 plasmid lp56'], ['Bacteria: Acholeplasma laidlawii PG-8A', 'Bacteria: Bacillus pseudofirmus OF4 plasmid pBpOF4-02'], ['Bacteria: Thermotoga maritima MSB8', 'Archaea: Pyrococcus furiosus DSM 3638'], ['Bacteria: Psychrobacter arcticus 273-4', 'Bacteria: Corynebacterium glutamicum ATCC 13032 DNA'], ['Bacteria: Thermobaculum terrenum ATCC BAA-798 chromosome 2', 'Bacteria: Deinococcus radiodurans R1 plasmid MP1'], ['Archaea: Sulfolobus solfataricus P2', 'Bacteria: Thermodesulfovibrio yellowstonii DSM 11347'], ['Archaea: Nitrosopumilus maritimus SCM1', 'Bacteria: Borrelia burgdorferi B31 plasmid lp56'], ['Bacteria: Alicyclobacillus acidocaldarius subsp. acidocaldarius DSM 446 plasmid pAACI03', 'Bacteria: Fibrobacter succinogenes subsp. succinogenes S85'], ['Bacteria: Pedobacter heparinus DSM 2366', 'Archaea: Methanococcus maripaludis C6'], ['Archaea: Methanosphaera stadtmanae DSM 3091', 'Archaea: Methanococcus aeolicus Nankai-3'], ['Bacteria: Fibrobacter succinogenes subsp. succinogenes S85', 'Archaea: Methanosarcina acetivorans str. C2A'], ['Bacteria: Gloeobacter violaceus PCC 7421 DNA', 'Archaea: Natronomonas pharaonis DSM 2160 plasmid PL23 complete genome'], ['Bacteria: Aquifex aeolicus VF5 plasmid ece1', 'Bacteria: Dictyoglomus thermophilum H-6-12'], ['Bacteria: Pirellula staleyi DSM 6068', 'Bacteria: Thermomicrobium roseum DSM 5159 plasmid'], ['Bacteria: Rhodopirellula baltica SH 1 complete genome', 'Bacteria: Chloroflexus aurantiacus J-10-fl'], ['Archaea: Haloterrigena turkmenica DSM 5511 plasmid pHTUR06', 'Archaea: Haloquadratum walsbyi DSM 16790 plasmid PL47 complete genome'], ['Archaea: Pyrococcus furiosus DSM 3638', 'Archaea: Pyrococcus horikoshii OT3 DNA'], ['Bacteria: Thermomicrobium roseum DSM 5159 plasmid', 'Bacteria: Bordetella pertussis strain Tohama I'], ['Archaea: Haloarcula marismortui ATCC 43049 plasmid pNG700', 'Archaea: Methanopyrus kandleri AV19'], ['Archaea: Hyperthermus butylicus DSM 5456', 'Archaea: Archaeoglobus fulgidus DSM 4304'], ['Bacteria: Magnetococcus marinus MC-1', 'Bacteria: Deinococcus radiodurans R1 plasmid MP1'], ['Archaea: Nanoarchaeum equitans Kin4-M', 'Bacteria: Dictyoglomus thermophilum H-6-12'], ['Archaea: Thermoplasma volcanium GSS1 DNA', 'Archaea: Methanococcus aeolicus Nankai-3'], ['Archaea: Methanococcoides burtonii DSM 6242', 'Bacteria: Rickettsia rickettsii str. Iowa'], ['Bacteria: Bacillus cereus Q1 plasmid pBc53', 'Bacteria: Bacillus pseudofirmus OF4 plasmid pBpOF4-02'], ['Bacteria: Elusimicrobium minutum Pei191', 'Bacteria: Chlorobium chlorochromatii CaD3'], ['Archaea: Methanocella paludicola SANAE DNA', 'Bacteria: Coxiella burnetii RSA 493 plasmid pQpH1'], ['Bacteria: Chloroflexus aurantiacus J-10-fl', 'Bacteria: Methylacidiphilum infernorum V4'], ['Archaea: Ignicoccus hospitalis KIN4/I', 'Archaea: Thermococcus gammatolerans EJ3'], ['Bacteria: Coxiella burnetii RSA 493 plasmid pQpH1', 'Bacteria: Bdellovibrio bacteriovorus complete genome'], ['Archaea: Halorhabdus utahensis DSM 12940', 'Archaea: Haloquadratum walsbyi DSM 16790 plasmid PL47 complete genome'], ['Archaea: Pyrococcus horikoshii OT3 DNA', 'Archaea: Pyrococcus furiosus DSM 3638'], ['Bacteria: Corynebacterium glutamicum ATCC 13032 DNA', 'Bacteria: Magnetococcus marinus MC-1'], ['Archaea: Caldivirga maquilingensis IC-167', 'Archaea: Metallosphaera sedula DSM 5348'], ['Bacteria: Clostridium acetobutylicum ATCC 824 megaplasmid pSOL1', 'Bacteria: Flavobacterium psychrophilum JIP02/86 complete genome'], ['Archaea: Thermofilum pendens Hrk 5 plasmid pTPEN01', 'Bacteria: Chlamydia trachomatis strain L2/434/Bu complete genome'], ['Bacteria: Akkermansia muciniphila ATCC BAA-835', 'Archaea: Cenarchaeum symbiosum A'], ['Bacteria: Deinococcus radiodurans R1 plasmid MP1', 'Bacteria: Bordetella pertussis strain Tohama I'], ['Bacteria: Thermanaerovibrio acidaminovorans DSM 6589', 'Bacteria: Thermus thermophilus HB8 plasmid pTT8 DNA'], ['Archaea: Sulfolobus acidocaldarius DSM 639', 'Archaea: Sulfolobus solfataricus P2'], ['Bacteria: Bordetella pertussis strain Tohama I', 'Bacteria: Deinococcus radiodurans R1 plasmid MP1'], ['Archaea: Halomicrobium mukohataei DSM 12286 plasmid pHmuk01', 'Bacteria: Rhodopirellula baltica SH 1 complete genome'], ['Bacteria: Thermodesulfovibrio yellowstonii DSM 11347', 'Archaea: Methanococcoides burtonii DSM 6242'], ['Bacteria: Prochlorococcus marinus str. AS9601', 'Archaea: Methanococcus aeolicus Nankai-3'], ['Bacteria: Bacteroides fragilis YCH46 plasmid pBFY46 DNA', 'Archaea: Methanocella paludicola SANAE DNA'], ['Bacteria: Rhizobium leguminosarum bv. trifolii WSM1325 plasmid pR132505', 'Bacteria: Pirellula staleyi DSM 6068'], ['Bacteria: Borrelia burgdorferi B31 plasmid lp56', 'Bacteria: Prochlorococcus marinus str. AS9601'], ['Archaea: Halorubrum lacusprofundi ATCC 49239 plasmid pHLAC01', 'Archaea: Haloquadratum walsbyi DSM 16790 plasmid PL47 complete genome'], ['Archaea: Aeropyrum pernix K1 DNA', 'Bacteria: Thermus thermophilus HB8 plasmid pTT8 DNA'], ['Archaea: Thermoplasma acidophilum DSM 1728 complete genome', 'Archaea: Sulfolobus solfataricus P2'], ['Archaea: Thermococcus gammatolerans EJ3', 'Archaea: Pyrococcus abyssi GE5 plasmid pGT5'], ['Archaea: Cenarchaeum symbiosum A', 'Bacteria: Acidobacterium capsulatum ATCC 51196'], ['Archaea: Methanococcus aeolicus Nankai-3', 'Archaea: Methanosphaera stadtmanae DSM 3091'], ['Bacteria: Thermus thermophilus HB8 plasmid pTT8 DNA', 'Bacteria: Thermanaerovibrio acidaminovorans DSM 6589'], ['Bacteria: Opitutus terrae PB90-1', 'Bacteria: Bordetella pertussis strain Tohama I'], ['Archaea: Methanosarcina acetivorans str. C2A', 'Bacteria: Chloroflexus aurantiacus J-10-fl'], ['Bacteria: Rhodospirillum rubrum ATCC 11170 plasmid', 'Archaea: Methanoculleus marisnigri JR1'], ['Bacteria: Bacillus pseudofirmus OF4 plasmid pBpOF4-02', 'Archaea: Metallosphaera sedula DSM 5348'], ['Bacteria: Synechococcus elongatus PCC 6301 DNA', 'Bacteria: Prochlorococcus marinus str. AS9601'], ['Archaea: Picrophilus torridus DSM 9790', 'Archaea: Methanococcus aeolicus Nankai-3'], ['Bacteria: Bdellovibrio bacteriovorus complete genome', 'Bacteria: Bacteroides fragilis YCH46 plasmid pBFY46 DNA'], ['Archaea: Staphylothermus marinus F1', 'Archaea: Thermoplasma volcanium GSS1 DNA'], ['Bacteria: Cyanothece sp. ATCC 51142 plasmid D', 'Archaea: Methanococcoides burtonii DSM 6242'], ['Bacteria: Carboxydothermus hydrogenoformans Z-2901', 'Bacteria: Mycoplasma genitalium G37 strain G-37 unique_12'], ['Bacteria: Gemmatimonas aurantiaca T-27 DNA', 'Bacteria: Opitutus terrae PB90-1'], ['Archaea: Pyrobaculum arsenaticum DSM 13514', 'Bacteria: Thermotoga maritima MSB8'], ['Archaea: Archaeoglobus profundus DSM 5631 plasmid pArcpr01', 'Archaea: Sulfolobus acidocaldarius DSM 639'], ['Archaea: Archaeoglobus fulgidus DSM 4304', 'Archaea: Methanosarcina acetivorans str. C2A'], ['Bacteria: Flavobacterium psychrophilum JIP02/86 complete genome', 'Archaea: Methanococcus maripaludis C6'], ['Archaea: Methanocaldococcus fervens AG86 plasmid pMEFER01', 'Archaea: Cenarchaeum symbiosum A'], ['Bacteria: Candidatus Amoebophilus asiaticus 5a2', 'Bacteria: Dictyoglomus thermophilum H-6-12'], ['Archaea: Methanobrevibacter ruminantium M1', 'Bacteria: Clostridium acetobutylicum ATCC 824 megaplasmid pSOL1'], ['Archaea: Pyrococcus abyssi GE5 plasmid pGT5', 'Archaea: Pyrococcus horikoshii OT3 DNA'], ['Archaea: Methanospirillum hungatei JF-1', 'Bacteria: Pedobacter heparinus DSM 2366'], ['Archaea: Methanococcus vannielii SB', 'Bacteria: Borrelia burgdorferi B31 plasmid lp56'], ['Archaea: Methanosaeta thermophila PT', 'Bacteria: Thermotoga maritima MSB8'], ['Bacteria: Campylobacter jejuni subsp. jejuni 81-176 plasmid pVir', 'Archaea: Methanococcus aeolicus Nankai-3']]

nodes = ['Archaea: Desulfurococcus kamchatkensis 1221n', 'Bacteria: Shewanella putrefaciens CN-32', 'Archaea: Methanosarcina mazei strain Goe1', 'Archaea: Natronomonas pharaonis DSM 2160 plasmid PL23 complete genome', 'Archaea: Methanocorpusculum labreanum Z', 'Bacteria: Methylacidiphilum infernorum V4', 'Bacteria: Acidobacterium capsulatum ATCC 51196', 'Archaea: Pyrobaculum aerophilum str. IM2', 'Archaea: Methanopyrus kandleri AV19', 'Archaea: Methanoculleus marisnigri JR1', 'Bacteria: Chlamydia trachomatis strain L2/434/Bu complete genome', 'Bacteria: Rickettsia rickettsii str. Iowa', 'Bacteria: Chlorobium chlorochromatii CaD3', 'Bacteria: Prochlorococcus marinus str. AS9601', 'Archaea: Haloquadratum walsbyi DSM 16790 plasmid PL47 complete genome', 'Archaea: Metallosphaera sedula DSM 5348', 'Bacteria: Dictyoglomus thermophilum H-6-12', 'Archaea: Picrophilus torridus DSM 9790', 'Bacteria: Mycoplasma genitalium G37 strain G-37 unique_12', 'Bacteria: Acholeplasma laidlawii PG-8A', 'Bacteria: Thermotoga maritima MSB8', 'Bacteria: Psychrobacter arcticus 273-4', 'Bacteria: Thermobaculum terrenum ATCC BAA-798 chromosome 2', 'Archaea: Sulfolobus solfataricus P2', 'Bacteria: Elusimicrobium minutum Pei191', 'Bacteria: Alicyclobacillus acidocaldarius subsp. acidocaldarius DSM 446 plasmid pAACI03', 'Bacteria: Pedobacter heparinus DSM 2366', 'Archaea: Methanosphaera stadtmanae DSM 3091', 'Bacteria: Fibrobacter succinogenes subsp. succinogenes S85', 'Archaea: Nitrosopumilus maritimus SCM1', 'Bacteria: Gloeobacter violaceus PCC 7421 DNA', 'Bacteria: Aquifex aeolicus VF5 plasmid ece1', 'Bacteria: Pirellula staleyi DSM 6068', 'Bacteria: Rhodopirellula baltica SH 1 complete genome', 'Archaea: Pyrococcus furiosus DSM 3638', 'Bacteria: Thermomicrobium roseum DSM 5159 plasmid', 'Archaea: Haloarcula marismortui ATCC 43049 plasmid pNG700', 'Archaea: Hyperthermus butylicus DSM 5456', 'Bacteria: Magnetococcus marinus MC-1', 'Archaea: Nanoarchaeum equitans Kin4-M', 'Archaea: Thermoplasma volcanium GSS1 DNA', 'Archaea: Methanococcoides burtonii DSM 6242', 'Bacteria: Bacillus cereus Q1 plasmid pBc53', 'Archaea: Methanocella paludicola SANAE DNA', 'Bacteria: Chloroflexus aurantiacus J-10-fl', 'Archaea: Ignicoccus hospitalis KIN4/I', 'Bacteria: Coxiella burnetii RSA 493 plasmid pQpH1', 'Archaea: Halorhabdus utahensis DSM 12940', 'Archaea: Pyrococcus horikoshii OT3 DNA', 'Bacteria: Corynebacterium glutamicum ATCC 13032 DNA', 'Archaea: Caldivirga maquilingensis IC-167', 'Bacteria: Clostridium acetobutylicum ATCC 824 megaplasmid pSOL1', 'Archaea: Thermofilum pendens Hrk 5 plasmid pTPEN01', 'Bacteria: Akkermansia muciniphila ATCC BAA-835', 'Bacteria: Deinococcus radiodurans R1 plasmid MP1', 'Bacteria: Thermanaerovibrio acidaminovorans DSM 6589', 'Archaea: Sulfolobus acidocaldarius DSM 639', 'Bacteria: Bordetella pertussis strain Tohama I', 'Archaea: Halomicrobium mukohataei DSM 12286 plasmid pHmuk01', 'Bacteria: Thermodesulfovibrio yellowstonii DSM 11347', 'Bacteria: Nostoc punctiforme PCC 73102 plasmid pNPUN05', 'Bacteria: Bacteroides fragilis YCH46 plasmid pBFY46 DNA', 'Bacteria: Rhizobium leguminosarum bv. trifolii WSM1325 plasmid pR132505', 'Bacteria: Borrelia burgdorferi B31 plasmid lp56', 'Archaea: Halorubrum lacusprofundi ATCC 49239 plasmid pHLAC01', 'Archaea: Aeropyrum pernix K1 DNA', 'Archaea: Thermoplasma acidophilum DSM 1728 complete genome', 'Archaea: Thermococcus gammatolerans EJ3', 'Archaea: Cenarchaeum symbiosum A', 'Archaea: Methanococcus aeolicus Nankai-3', 'Bacteria: Thermus thermophilus HB8 plasmid pTT8 DNA', 'Bacteria: Carboxydothermus hydrogenoformans Z-2901', 'Bacteria: Opitutus terrae PB90-1', 'Bacteria: Bacillus pseudofirmus OF4 plasmid pBpOF4-02', 'Bacteria: Rhodospirillum rubrum ATCC 11170 plasmid', 'Archaea: Methanosarcina acetivorans str. C2A', 'Bacteria: Synechococcus elongatus PCC 6301 DNA', 'Archaea: Methanococcus maripaludis C6', 'Bacteria: Bdellovibrio bacteriovorus complete genome', 'Archaea: Staphylothermus marinus F1', 'Bacteria: Cyanothece sp. ATCC 51142 plasmid D', 'Archaea: Haloterrigena turkmenica DSM 5511 plasmid pHTUR06', 'Bacteria: Gemmatimonas aurantiaca T-27 DNA', 'Archaea: Pyrobaculum arsenaticum DSM 13514', 'Archaea: Archaeoglobus profundus DSM 5631 plasmid pArcpr01', 'Archaea: Archaeoglobus fulgidus DSM 4304', 'Bacteria: Flavobacterium psychrophilum JIP02/86 complete genome', 'Archaea: Methanocaldococcus fervens AG86 plasmid pMEFER01', 'Bacteria: Candidatus Amoebophilus asiaticus 5a2', 'Archaea: Methanobrevibacter ruminantium M1', 'Archaea: Pyrococcus abyssi GE5 plasmid pGT5', 'Archaea: Methanospirillum hungatei JF-1', 'Archaea: Methanococcus vannielii SB', 'Archaea: Methanosaeta thermophila PT', 'Bacteria: Campylobacter jejuni subsp. jejuni 81-176 plasmid pVir']

types = ['Archaea', 'Bacteria', 'Archaea', 'Archaea', 'Archaea', 'Bacteria', 'Bacteria', 'Archaea', 'Archaea', 'Archaea', 'Bacteria', 'Bacteria', 'Bacteria', 'Bacteria', 'Archaea', 'Archaea', 'Bacteria', 'Archaea', 'Bacteria', 'Bacteria', 'Bacteria', 'Bacteria', 'Bacteria', 'Archaea', 'Archaea', 'Bacteria', 'Bacteria', 'Archaea', 'Bacteria', 'Bacteria', 'Bacteria', 'Bacteria', 'Bacteria', 'Archaea', 'Archaea', 'Bacteria', 'Archaea', 'Archaea', 'Bacteria', 'Archaea', 'Archaea', 'Archaea', 'Bacteria', 'Bacteria', 'Archaea', 'Bacteria', 'Archaea', 'Bacteria', 'Archaea', 'Archaea', 'Bacteria', 'Archaea', 'Bacteria', 'Archaea', 'Bacteria', 'Bacteria', 'Bacteria', 'Archaea', 'Bacteria', 'Archaea', 'Bacteria', 'Bacteria', 'Bacteria', 'Bacteria', 'Bacteria', 'Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea', 'Bacteria', 'Bacteria', 'Archaea', 'Bacteria', 'Bacteria', 'Bacteria', 'Archaea', 'Bacteria', 'Archaea', 'Bacteria', 'Bacteria', 'Bacteria', 'Archaea', 'Archaea', 'Archaea', 'Bacteria', 'Archaea', 'Bacteria', 'Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea', 'Bacteria']

reseau2(relations, nodes, types)        
"""
path = "D:\\Documentos\\Master\\Deuxieme_annee\\GENOM\\Projet\\Projet 2\\genomes"
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
message = "Liste des genomes charges: \n"
sequences  = {}
cont = 1
for f in onlyfiles:
    s = path+'/'+f
    name, genome = lireSeq(s)
    print (cont)
    cont+=1
"""    
"""    
reseau2(relations, Nodes, c)        
relations = [['Archaea: Nanoarchaeum equitans Kin4-M', 'Archaea: Methanosarcina mazei strain Goe1'], ['Archaea: Methanosarcina mazei strain Goe1', 'Archaea: Nanoarchaeum equitans Kin4-M'], ['Archaea: Archaeoglobus fulgidus DSM 4304', 'Archaea: Methanosarcina mazei strain Goe1']]
Nodes = ['Archaea: Nanoarchaeum equitans Kin4-M', 'Archaea: Methanosarcina mazei strain Goe1', 'Archaea: Archaeoglobus fulgidus DSM 4304']    
c = ['Archaea', 'Archaea', 'Archaea']
from igraph import *
g = Graph(directed=True)
g.add_vertices(Nodes)
g.vs['label'] = Nodes

for i in relations:
    g.add_edges([(i[0],i[1])])

layout = g.layout_kamada_kawai()
color_dict = {"Archaea": "blue", "Bacteria": "pink"}
g.vs["color"] = [color_dict[t] for t in c]
plot(g, layout)#, bbox = (300, 300), margin = 20)
""" 

"""    
    print (matrice)
    f=open("Matrice.txt",'w')
    for i in range(len(matrice)):
        for j in range(len(matrice)):
            f.write(matrice[i,j]+'\t') 
    f.close()
    """
    #distance_matrix = r.matrice()
        
    #genome_hclust = r.hclust(distance_matrix)
    #r.plot(genome_hclust)
"""
m = np.array([[  0,   1.02047102e-04,   7.04264773e-05],
 [  1.02047102e-04,   0,   5.44397724e-05],
 [  7.04264773e-05,   5.44397724e-05, 0]])

f=open("Matrice.txt",'w')
for i in range(len(m)):
    for j in range(len(m)):
        print (m[i,j])
        f.write(str(m[i,j])+' ') 
    f.write('\n') 
f.close()
"""

"""
#print (path_dir)
k=6
name, seq = lireSeq(path_dir+"\\seq.fasta")
alphabet_nt={'A':0,'T':0,'G':0,'C':0}

dico_k_mers = do_dic_kmers(k, alphabet_nt)
kmers = find_kmers(dico_k_mers, k, seq)
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
