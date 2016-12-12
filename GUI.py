"""
MAIN - INTERFACE GRAPHIQUE
"""
import random
import sys, os
from  background import *
# choix de la version de python
version = sys.version_info[0]
from os import listdir
from os.path import isfile, join
path_dir = os.getcwd()
import shutil
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_pdf import PdfPages
import timeit
import subprocess

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 8}

plt.rc('font', **font)


if os.path.exists(path_dir+ "/images" ):
    shutil.rmtree(path_dir+ "/images" , ignore_errors=True)
    os.makedirs(path_dir+ "/images")
else:
    os.makedirs(path_dir+ "/images")

if version == 2:
	#### pour python 2
     import Tkinter as tkinter
     from Tkinter import Tk, Text, BOTH, W, N, E, S, END, CENTER
     from ttk import Frame, Button, Label, Style
     import tkFileDialog as filedialog
     from ScrolledText import ScrolledText
     from tkFont import Font
     from tkMessageBox import Message
     import tkMessageBox as mbox

if version == 3:
	#### pour python 3
    import tkinter as tkinter
    from tkinter import *
    from tkinter import Tk, Text, BOTH, W, N, E, S, END, CENTER
    from tkinter.ttk import Frame, Button, Label, Style
    import tkinter.filedialog as filedialog
    from tkinter.scrolledtext import ScrolledText  
    from tkinter.font import Font
    from tkinter.messagebox import Message
    import tkinter.messagebox as mbox
    import tkinter.ttk as ttk


alphabet_nt={'A':0,'T':0,'G':0,'C':0}
 
class Kmers(Frame):
    
    def __init__(self, parent):
        Frame.__init__(self, parent)   
        
        self.customFont = Font(family="consolas", size=10)
        self.seq = ""
        self.kmer = 2
        self.parent = parent
        self.initUI()
        
        
    def initUI(self):
      
        self.parent.title("Kmers")
        self.pack(fill=BOTH, expand=False)
        
        Label(self, text="Fichier de la Sequence:").grid(sticky=W+E+S+N, row=0, column=0, padx =5, pady=5)
        Label(self, text="Nombre des Kmers:").grid(sticky=W+E+S+N, row=1, column=0, padx =5, pady=5)

        self.btn = Button(self, text="Open File", command = self.onOpen)
        self.btn.grid(sticky=E, row=0, column=0, padx=5, pady=5)
        
        self.e1 = tkinter.Entry(self, state="disabled")
        self.e1.configure(width=12)
        self.e1.grid(sticky=E, row=1, column=0, padx =5, pady=5)
        
        self.btn2 = Button(self, text="Calcule", state="disabled", command = self.onCalcule)
        self.btn2.grid(sticky=W+E+S+N, row=1, column=1, padx=5, pady=5)
        
        lbl = Label(self, text="Infos")
        #lbl.grid(sticky=W, row=2, column=0, padx=5)
        
        self.v = tkinter.IntVar()
        self.v = 2
        r1 = tkinter.Radiobutton(self, text="Euclidenne", variable=self.v, value=1)
        r2 = tkinter.Radiobutton(self, text="Correlation", variable=self.v, value=2)
        
        r1.grid(sticky=W, row=2, column=0, padx=5)
        r2.grid(sticky=W, row=2, column=1, padx=5)
        
        lbl2 = Label(self, text="List Kmers")
        lbl2.grid(sticky=W, row=0, column=2, padx=5)
        
        self.area = ScrolledText(self, wrap=tkinter.WORD)
        self.area.configure(state="disabled", font=self.customFont, height= 20, width=50)
        self.area.grid(padx =5, pady=2, row=3, columnspan = 2, rowspan=30)
        
        self.listbox = tkinter.Listbox(self)
        self.listbox.config(height= 25, width=58)
        self.listbox.grid(pady=5, row=1, column = 2,  rowspan=30)
        
        lbl = Label(self, text="1. Proportion genome")
        lbl.grid(sticky=W+E+S+N, row=0, column=4, padx=5)
        
        lbl = Label(self, text="Pas:")
        lbl.grid(sticky=W+E+S+N, row=1, column=4, padx=5)
        
        lbl = Label(self, text="Fenetre:")
        lbl.grid(sticky=W+E+S+N, row=2, column=4, padx=5)
        
        self.e2 = tkinter.Entry(self, state="disabled")
        self.e2.configure(width=12)
        self.e2.grid(sticky=E, row=1, column=4, padx =5, pady=5)
        
        self.e3 = tkinter.Entry(self, state="disabled")
        self.e3.configure(width=12)
        self.e3.grid(sticky=E, row=2, column=4, padx =5, pady=5)
        
        self.btn3 = Button(self, text="Calcule", state="disabled", command = self.onCalcule2)
        self.btn3.grid(sticky=W+E+S+N, row=2, column=5, padx=5, pady=5)
        
        lbl = Label(self, text="Calcule subparties")
        lbl.grid(sticky=W+E+S+N, row=3, column=4, padx=5)
        
        lbl = Label(self, text="Pas:")
        lbl.grid(sticky=W+E+S+N, row=4, column=4, padx=5)
        
        lbl = Label(self, text="Fenetre:")
        lbl.grid(sticky=W+E+S+N, row=5, column=4, padx=5)
        
        self.e4 = tkinter.Entry(self, state="disabled")
        self.e4.configure(width=12)
        self.e4.grid(sticky=E, row=4, column=4, padx =5, pady=5)
        
        self.e5 = tkinter.Entry(self, state="disabled")
        self.e5.configure(width=12)
        self.e5.grid(sticky=E, row=5, column=4, padx =5, pady=5)
        
        self.btn4 = Button(self, text="Calcule", state="disabled", command = self.onCalcule3)
        self.btn4.grid(sticky=W+E+S+N, row=5, column=5, padx=5, pady=5)
        
        self.imp = tkinter.BooleanVar()    
        self.cb = tkinter.Checkbutton(self, text="Imprimer Partie Genome diferente", state="disabled", variable=self.imp, command=self.onClick)
        self.cb.grid(sticky=W+E+S+N, row=6, column=4, columnspan = 2, padx=5, pady=5)
        
        self.btn5 = Button(self, text="Creer Phylo", state="disabled", command = self.onCalcule4)
        self.btn5.grid(sticky=W+E+S+N, row=7, column=4, columnspan = 2, padx=5, pady=5)
        
        self.btn6 = Button(self, text="Creer Reseau Relations", state="disabled", command = self.onCalcule5)
        self.btn6.grid(sticky=W+E+S+N, row=8, column=4, columnspan = 2, padx=5, pady=5)
        
        
        lbl = Label(self, text="2. Choisir genome / genome Predicte")
        lbl.grid(sticky=W+E+S+N, row=0, column=7, columnspan = 2, padx=5)
        
        lbl = Label(self, text="Taille morceau")
        lbl.grid(sticky=W+E+S+N, row=1, column=7, padx=5)
        
        self.e7 = tkinter.Entry(self, state="disabled")
        self.e7.configure(width=12)
        self.e7.grid(sticky=E, row=1, column=8, padx =5, pady=5)
        
        self.btn7 = Button(self, text="Cherching morceau", state="disabled", command = self.onCalcule6)
        self.btn7.grid(sticky=W+E+S+N, row=2, column=7, columnspan = 2, padx=5, pady=5)
        
        lbl = Label(self, text="3. Experience aleatoire un genome")
        lbl.grid(sticky=W+E+S+N, row=4, column=7, columnspan = 2, padx=5)
                
        lbl = Label(self, text="Taille morceau")
        lbl.grid(sticky=W+E+S+N, row=5, column=7, padx=5)
        
        self.e8 = tkinter.Entry(self, state="disabled")
        self.e8.configure(width=12)
        self.e8.grid(sticky=E, row=5, column=8, padx =5, pady=5)
        
        lbl = Label(self, text="Nombre de fois")
        lbl.grid(sticky=W+E+S+N, row=6, column=7, padx=5)
        
        self.e9 = tkinter.Entry(self, state="disabled")
        self.e9.configure(width=12)
        self.e9.grid(sticky=E, row=6, column=8, padx =5, pady=5)
        
        self.btn8 = Button(self, text="Lancer l'experience ", state="disabled", command = self.onCalcule7)
        self.btn8.grid(sticky=W+E+S+N, row=7, column=7, columnspan = 2, padx=5, pady=5)
        
        
        lbl = Label(self, text="3. Experience aleatoire tous les genomes")
        lbl.grid(sticky=W+E+S+N, row=8, column=7, columnspan = 2, padx=5)
                
        lbl = Label(self, text="Taille morceau")
        lbl.grid(sticky=W+E+S+N, row=9, column=7, padx=5)
        
        self.e10 = tkinter.Entry(self, state="disabled")
        self.e10.configure(width=12)
        self.e10.grid(sticky=E, row=9, column=8, padx =5, pady=5)
        
        lbl = Label(self, text="Nombre de fois")
        lbl.grid(sticky=W+E+S+N, row=10, column=7, padx=5)
        
        self.e11 = tkinter.Entry(self, state="disabled")
        self.e11.configure(width=12)
        self.e11.grid(sticky=E, row=10, column=8, padx =5, pady=5)
        
        self.btn9 = Button(self, text="Lancer l'experience ", state="disabled", command = self.onCalcule8)
        self.btn9.grid(sticky=W+E+S+N, row=11, column=7, columnspan = 2, padx=5, pady=5)
        
        
    def onOpen(self):
        ftypes = [('Fna files', '*.fna')]
        dlg = filedialog.Open (self, filetypes = ftypes)
        fl = dlg.show()
        n = len(fl.split('/')[-1])+1
        path = fl[:-n]
        start = timeit.default_timer()

        onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
        message = "Liste des genomes charges: \n"
        self.sequences  = {}
        for f in onlyfiles:
            s = path+'/'+f
            name, genome = lireSeq(s)
            self.sequences[name] = genome
            message = message + name+ '\n'
        stop = timeit.default_timer()
        time_total = stop-start

        message += '\n'
        message += ' --> Time reading : ' + str(time_total) +'\n\n'
        self.area.configure(state="normal")
        self.area.insert(END, message)
        
        self.area.configure(state="disabled") 
        self.btn2.configure(state="normal")      
        self.e1.configure(state="normal")

    def onCalcule(self):
        k = self.e1.get()
        try:
            self.kmer = int(k)
            self.e1.configure(state="disabled")
            self.btn2.configure(state="disabled")
            
            self.area.configure(state="normal")
            message = "Kmer calcule : "+ str(k)+"\n"
            self.area.insert(END, message)
            self.area.configure(state="disabled") 
            
            start = timeit.default_timer()
    
            self.dico_k_mers = do_dic_kmers(self.kmer, alphabet_nt)
            try:
                self.dicProfils = do_profils_all_genomes(self.sequences, self.dico_k_mers, self.kmer)
                """
                cont = 0
                for key in self.dicProfils.keys():
                    l = list(self.dicProfils[key].values()).count(0)
                    string = key+ "  Lenght Profil: "+ str(len(self.dicProfils[key])-l)
                    self.listbox.insert(cont, string)
                    cont += 1
                """
            except Exception as e:
                print ("ici INI",e)
            
            stop=timeit.default_timer()
            time_total=stop-start

            message = ' --> Time calculating Profils : ' + str(time_total) +'\n\n'
            
            self.area.configure(state="normal")
            self.area.insert(END, message)
            self.area.configure(state="disabled") 
            
            self.btn3.configure(state="normal")      
            self.e2.configure(state="normal")
            self.e3.configure(state="normal")
            
            self.btn7.configure(state="normal")      
            self.e7.configure(state="normal")
            self.btn8.configure(state="normal")      
            self.e8.configure(state="normal")
            self.e9.configure(state="normal")
            self.btn9.configure(state="normal")      
            self.e10.configure(state="normal")
            self.e11.configure(state="normal")
            
        
        except Exception as e:
            print ("ici INI2",e)
            mbox.showwarning("Warning", "Number must be Integer")
    
    def onCalcule2(self):
        pas = self.e2.get()
        fenetre = self.e3.get()
        self.pas = []
        self.fenetres = []
        self.positions = {}
        self.infoDistances = {}
        self.partieGenomeDif = {}
        self.cont = 0

        try:
            self.cont = 0
            self.pas.append(int(pas))
            self.fenetres.append(int(fenetre))
            
            self.area.configure(state="normal")
            message = "Proportion Genome : pas: "+ str(pas)+" fenetre: "+str(fenetre)+" \n"
            message += "Calcule distance..."+"\n"
            self.area.insert(END, message)
            self.area.configure(state="disabled") 
            
            start = timeit.default_timer()
            
            self.parties = do_sub_parts(self.dicProfils, self.sequences, int(pas), int(fenetre), self.dico_k_mers, self.kmer, self.v)
            
            pdfFile = path_dir+'/images/FirstParties.pdf'
            pp = PdfPages(pdfFile)
            
            means = do_mean_sub_parts(self.parties)
            ecarts = do_sd_sub_parts(self.parties, means)
            
            message = "Position dans le genome avec la max dist: \n"
            try:
                for i in self.parties.keys():
                    self.positions[i] = [self.parties[i].index(max(self.parties[i]))*(int(pas))]
                    listI = []
                    listI.append([max(self.parties[i]), means[i]+ecarts[i]*2])
                    self.infoDistances[i] = listI
                    message += i+' pos: '+ str(self.positions[i][self.cont])+"\n"
                message += '\n'
            except Exception as e:
                print ("ici 2",e)
            
            
            try:
                self.plot_for_each_genome(self.parties, pp, self.pas, self.positions, self.fenetres, self.cont, True , means, ecarts)
            except Exception as e:
                print ("ici 1",e)
            
            subprocess.Popen([pdfFile],shell=True)
            
            self.area.configure(state="normal")
            self.area.insert(END, message)
            self.area.configure(state="disabled")
            
            try:
                self.partieGenomeDif = get_partieDif_all_genomes(self.sequences, self.positions, self.pas, int(fenetre), self.cont)
            except Exception as e:
                print ("ici 3",e)
            
            stop=timeit.default_timer()
            time_total=stop-start

            message = ' --> Time calculating sub-parties : ' + str(time_total) +'\n\n'
            
            self.area.configure(state="normal")
            self.area.insert(END, message)
            self.area.configure(state="disabled") 
            
            self.btn3.configure(state="disabled")      
            self.e2.configure(state="disabled")
            self.e3.configure(state="disabled")
            
            self.btn4.configure(state="normal")      
            self.e4.configure(state="normal")
            self.e5.configure(state="normal")
            
        except Exception as e:
            print ("ici 3",e)
            mbox.showwarning("Warning", "Number must be Integer")                
    
            
    def myplotcode(self, x, y):
        fig = plt.figure(figsize=(40, 40))
        ax = fig.add_subplot(111)
        x_points = range(0,len(x))
        y_points = y
        ax.plot(x_points, y_points, 'b')
        ax.set_xlabel('Parties Genome')
        ax.set_ylabel('Distance')
        ax.set_title('Distribution proportion Genome')
        return fig
    
    def plot_for_each_genome(self, sub_parts_for_each_genome, pp, pas, positions, fenetres, cont, boolean, means, ecarts):
        try:
            for i in sub_parts_for_each_genome.keys():
                fig = plt.figure()
                ax = fig.add_subplot(111)
                pos = 0
                position = positions[i]
                if boolean:   
                    pos = int(position[cont]-(fenetres[cont]/2))
                else:
                    pos = int(position[-1]-(fenetres[-2]/2))
                    if pos < 0:
                        pos = 0
                if cont == 0:
                    xVector = range(0,len(sub_parts_for_each_genome[i])*pas[cont], pas[cont])
                    x_points = xVector
                else:
                    xVector = range(pos, pos + len(sub_parts_for_each_genome[i])*pas[cont], pas[cont])
                    x_points = xVector
                
                y_points = sub_parts_for_each_genome[i]
                plt.plot(x_points, y_points, 'b')
                ax.set_xlabel('Position')
                ax.set_ylabel('Distance')
                ax.set_title(i)
                x = x_points
                y = len(x)*[means[i]]
                plt.plot(x, y)
                x = x_points
                y = len(x)*[means[i]+ecarts[i]*2]
                plt.plot(x, y)
                x = x_points
                y = len(x)*[means[i]-ecarts[i]*2]
                plt.plot(x, y)
                pp.savefig()
            pp.close()
        except Exception as e:
                print (e)
    

    def onCalcule3(self):
        pas = self.e4.get()
        fenetre = self.e5.get()
        self.cont += 1
        try:
            if int(fenetre) == self.fenetres[-1] or int(fenetre) > self.fenetres[-1]:
                mbox.showwarning("Warning", "La fenetre foit etre de taille inferieure a la derniere fenetre")
                return
            else:
                self.pas.append(int(pas))
                self.fenetres.append(int(fenetre))
                
            self.area.configure(state="normal")
            message = "Proportion subpartie : pas: "+ pas+" fenetre: "+ fenetre+" \n"
            message += "Calcule distance..."+"\n"
            self.area.insert(END, message)
            self.area.configure(state="disabled") 
            
            start = timeit.default_timer()
            
            self.parties = do_sub_parts(self.dicProfils, self.partieGenomeDif, int(pas), int(fenetre), self.dico_k_mers, self.kmer, self.v)
            
            means = do_mean_sub_parts(self.parties)
            ecarts = do_sd_sub_parts(self.parties, means)
            
            pdfFile = path_dir+'/images/SubParties'+str(self.cont)+'.pdf'
            pp = PdfPages(pdfFile)
            try:
                message = "Position dans le genome avec la max dist: \n"

                for i in self.parties.keys():
                    list_i = self.positions[i]
                    pos = 0
                    pos = int((list_i[-1] - (self.fenetres[-2])/2) + self.parties[i].index(max(self.parties[i]))*int(pas))
                    if pos < 0 :
                        pos = 0
                    list_i.append(pos)
                    self.positions[i] = list_i
                    listj = self.infoDistances[i]
                    listj.append([max(self.parties[i]), means[i]+ecarts[i]*2])
                    self.infoDistances[i] = listj
                    message += i+' pos: '+ str(pos)+'\n'
                message += '\n'
            except Exception as e:
                print (e)
            
            try:
                self.plot_for_each_genome(self.parties, pp, self.pas, self.positions, self.fenetres, self.cont, False, means, ecarts)
            except Exception as e:
                print ("ici 1",e)
            subprocess.Popen([pdfFile],shell=True)    
            
            self.area.configure(state="normal")
            self.area.insert(END, message)
            self.area.configure(state="disabled")
            
            try:
                self.partieGenomeDif = get_partieDif_all_genomes(self.sequences, self.positions, self.pas, int(fenetre), self.cont)
                matrice, self.listCles = transfert_study(self.dicProfils, self.partieGenomeDif, self.dico_k_mers, self.kmer, self.v)    
                liste_probable_transfert = tranfert_study_result(matrice, self.listCles, self.infoDistances)
                self.relations = file_result(liste_probable_transfert, "Resultats.txt")
                self.types = []
                for i in self.listCles:
                    self.types.append(i.split(":")[0])
            except Exception as e:
                print ("ici 4",e)
            
            stop=timeit.default_timer()
            time_total=stop-start

            message = ' --> Time calculating sub-parties : ' + str(time_total) +'\n\n'
            
            self.area.configure(state="normal")
            self.area.insert(END, message)
            self.area.configure(state="disabled") 
            
            self.cb.configure(state="normal")
            self.btn5.configure(state="normal")      
            self.btn6.configure(state="normal")  
        except Exception as e:
            print ("ici 4",e)
            mbox.showwarning("Warning", "Number must be Integer")                
    
    def onClick(self):
        if self.imp.get() == True:
            self.area.configure(state="normal")
            self.area.insert(END, "Fichier creee et sauvegarde dans le repertoire actuel : Nom fichier : SeqDifs.txt \n\n")
            self.area.configure(state="disabled")   
            sequencesDif(self.partieGenomeDif)
    
    def onCalcule4(self):
        pp = PdfPages(path_dir+'\\images\\Phylogenie.pdf')
        matrice, liste = do_data_matrix(self.dicProfils)
        Z = linkage(matrice, 'ward')
        fig = plt.figure(figsize=(30, 10))
        ax = fig.add_subplot(111)
        ax.set_title('Hierarchical Clustering Dendrogram Genomes')  
        ax.set_xlabel('Distance')
        ax.set_ylabel('Genomes')
        dendrogram(Z, leaf_rotation=0., leaf_font_size=6., orientation='right', labels = liste)
        fig.show()
        pp.savefig()
        pp.close()
        self.fenetrePlot(fig, "Phylogenie" )
        
    def fenetrePlot(self, fig, title):
        self.plot = tkinter.Toplevel(self)
        self.plot.geometry("700x700+100+100")
        self.plot.wm_title(title)
        self.plot.canvas = FigureCanvasTkAgg(fig, master=self.plot)
        self.plot.canvas.show()
        self.plot.canvas.get_tk_widget().pack()
    
    def onCalcule5(self):
        reseau2(self.relations, self.listCles, self.types)
        #self.fenetrePlot(graph, title)
        try: 
            subprocess.Popen([path_dir+"\\Reseau.pdf"],shell=True)
        except Exception as e:
            print ("ici 4",e) 
    
    def onCalcule6(self):
        taille = self.e7.get()
        try: 
            self.taille = int (taille)
            
            key = random.choice(list(self.sequences.keys()))
            genome = self.sequences[key]
            rand = random.randint(0, len(genome) - self.taille)
            randomSeq = genome[rand :  rand+self.taille]
            message = 'Genome Aleatoire : \n'
            message += key+'\n'
            
            self.area.configure(state="normal")
            self.area.insert(END, message)
            self.area.configure(state="disabled")
            
            profilMorceau = do_profil_genome(randomSeq, self.dico_k_mers, self.kmer)
            ind = transfert_studyBoutGenome(self.dicProfils, profilMorceau, self.v)
            
            message = 'Genome predicte \n'
            message += ind +'\n\n'
            self.area.configure(state="normal")
            self.area.insert(END, message)
            self.area.configure(state="disabled")
            
        except Exception as e:
            print ("ici 4",e) 
            mbox.showwarning("Warning", "Number must be Integer")
    
    def onCalcule7(self):
        taille = self.e8.get()
        expFois = self.e9.get()
        try:
            self.taille = int(taille)
            self.fois = int(expFois)            
            
            vectResultat, genome = experiment(self.dicProfils, self.sequences, self.dico_k_mers, self.kmer, self.fois, self.taille, self.v)
            erreur = errorRate(vectResultat)
            
            message = 'Genome Aleatoire\n'
            message += genome +'\n'
            message += "Taux d'erreur pour le genome choisi : " + str(erreur) +'\n\n'
            self.area.configure(state="normal")
            self.area.insert(END, message)
            self.area.configure(state="disabled")
        except Exception as e:
            print ("ici 4",e) 
            mbox.showwarning("Warning", "Number must be Integer")
    
    def onCalcule8(self):
        taille = self.e10.get()
        expFois = self.e11.get()
        try:
            self.taille = int(taille)
            self.fois = int(expFois)            
            
            erreur, genomes = experiment2(self.dicProfils, self.sequences, self.dico_k_mers, self.kmer, self.fois, self.taille, self.v)
            types = []
            for i in genomes:
                types.append(i.split(":")[0])
            fig = plt.figure(figsize=(25, 10))
            width = 0.2
            ax = fig.add_subplot(111)
            ax.set_title("Taux d'erreur predictions")  
            ax.set_xlabel('Genomes')
            ax.set_ylabel("Taux d'erreur")
            ind = np.arange(len(genomes)) 
            plt.bar(ind, erreur, width, color='b')
            plt.xticks(ind + width/2., types, rotation=45)
            plt.yticks(np.arange(0, 1.5, 0.2))
            fig.show()
            self.fenetrePlot(fig, "Taux d'erreur" )
            
            message = 'Genome aleatoire\n'
            message += 'Calcule pour tous les genomes\n'
            message += 'Prediction des '+str(self.fois)+" bouts aleatoires pour chaque genome \n\n"
            self.area.configure(state="normal")
            self.area.insert(END, message)
            self.area.configure(state="disabled")
        except Exception as e:
            print ("ici 4",e) 
            mbox.showwarning("Warning", "Number must be Integer")
            
    
        
def main():
    root = Tk()
    root.geometry("1200x500+200+200")
    app = Kmers(root)
    root.mainloop()  


if __name__ ==  "__main__" :
    main()
