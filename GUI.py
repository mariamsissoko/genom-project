"""
MAIN - INTERFACE GRAPHIQUE
"""

import sys, os
from mod_background import *
# choix de la version de python
version = sys.version_info[0]

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
     import numpy as np
     from matplotlib.figure import Figure
     import matplotlib.pyplot as plt
     from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


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

class DockingProteines(Frame):
    
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
        lbl.grid(sticky=W, row=2, column=0, padx=5)
        
        lbl2 = Label(self, text="List Kmers")
        lbl2.grid(sticky=W, row=0, column=2, padx=5)
        
        self.area = ScrolledText(self, wrap=tkinter.WORD)
        self.area.configure(state="disabled", font=self.customFont, height= 20, width=50)
        self.area.grid(padx =5, pady=2, row=3, columnspan = 2, rowspan=30)
        
        self.listbox = tkinter.Listbox(self)
        self.listbox.config(height= 25, width=40)
        self.listbox.grid(pady=5, row=1, column = 2, rowspan=30)
        
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
        
        lbl = Label(self, text="Calcule de sub-parties")
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
        self.cb = tkinter.Checkbutton(self, text="Imprimer Partie Genome dif", state="disabled", variable=self.imp, command=self.onClick)
        self.cb.grid(sticky=W+E+S+N, row=6, column=4, padx=5, pady=5)
        
    def onOpen(self):
        ftypes = [('Fasta files', '*.fasta')]
        dlg = filedialog.Open (self, filetypes = ftypes)
        fl = dlg.show()
        self.seq = lireSeq(fl)
        self.area.configure(state="normal")
        message = "Lenght sequence : "+ str(len(self.seq))+"\n"
        self.area.insert(END, message)
        self.area.configure(state="disabled") 
        self.btn2.configure(state="normal")      
        self.e1.configure(state="normal")
        self.btn.configure(state="disabled") 
        self.btn3.configure(state="normal")      
        self.e2.configure(state="normal")
        self.e3.configure(state="normal")
        

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
            self.dico_k_mers = do_dic_kmers(self.kmer, alphabet_nt)
            self.kmers = find_kmers(self.dico_k_mers, self.kmer, self.seq)
            self.propotionGenome = proportion_of_kmers(len(self.seq)-self.kmer+1, self.kmers)
            cont = 0
            for key in self.kmers:
                #"string = '{} {:>10}{:>10}'.format( key, str(self.kmers[key]), str(self.propotionGenome[key]))
                string = key+ "   -   "+ str(self.kmers[key]) +"   -   "+ str(round(self.propotionGenome[key],6))
                self.listbox.insert(cont, string)
                cont += 1

        except:
            mbox.showwarning("Warning", "Number must be Integer")
    
    def onCalcule2(self):
        pas = self.e2.get()
        fenetre = self.e3.get()
        self.pas = []
        self.positions = []
        try:
            self.cont = 0
            self.pas.append(int(pas))
            self.fenetre = int(fenetre)
            
            self.cb.configure(state="normal")
            
            self.area.configure(state="normal")
            message = "Proportion Genome : pas: "+ str(pas)+" fenetre: "+str(self.fenetre)+" \n"
            self.area.insert(END, message)
            self.area.configure(state="disabled") 
            
            self.parties = proportion_along_genome(self.seq, int(pas), self.fenetre, self.dico_k_mers, self.kmer) 
            self.area.configure(state="normal")
            message = "Nombre parties selectionées : "+ str(len(self.parties))+" \n"
            self.area.insert(END, message)
            self.area.configure(state="disabled")
            
            self.distances = calculeDistanceParties(self.propotionGenome, self.parties)
            self.position = self.distances.index(max(self.distances))
            self.positions.append(self.position)
            
            self.area.configure(state="normal")
            message = "Calcule distance Euclidienne"+"\n"
            self.area.insert(END, message)
            self.area.configure(state="disabled")
            
            self.area.configure(state="normal")
            message = "Position dans le genome avec la max dist: "+ str(self.position)+"\n\n"
            self.area.insert(END, message)
            self.area.configure(state="disabled")
            
            self.partieGenomeDif = self.seq[self.position * int(pas) : self.position * int(pas) + self.fenetre]
            
            self.btn3.configure(state="disabled")      
            self.e2.configure(state="disabled")
            self.e3.configure(state="disabled")
            
            self.btn4.configure(state="normal")      
            self.e4.configure(state="normal")
            self.e5.configure(state="normal")
                
            self.fig = self.myplotcode(self.parties, self.distances)
            self.fenetrePlot()

        except:
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
    
    def fenetrePlot(self):
        self.plot = tkinter.Toplevel(self)
        self.plot.geometry("500x500+300+200")
        self.plot.wm_title("Plot distribution proportions Genome")
        self.plot.canvas = FigureCanvasTkAgg(self.fig, master=self.plot)
        self.plot.canvas.show()
        self.plot.canvas.get_tk_widget().pack()

    def onCalcule3(self):
        pas = self.e4.get()
        fenetre = self.e5.get()
        try:
            self.pas.append(int(pas))
            self.fenetre2 = int(fenetre)
            
            self.area.configure(state="normal")
            message = "Proportion subpartie : pas: "+ str(pas)+" fenetre: "+str(self.fenetre2)+" \n"
            self.area.insert(END, message)
            self.area.configure(state="disabled") 
            
            self.subparties = proportion_along_genome(self.partieGenomeDif, int(pas), self.fenetre2, self.dico_k_mers, self.kmer) 
            
            self.area.configure(state="normal")
            message = "Nombre parties selectionées : "+ str(len(self.subparties))+" \n"
            self.area.insert(END, message)
            self.area.configure(state="disabled")
            
            self.distances2 = calculeDistanceParties(self.propotionGenome, self.subparties)
            self.position2 = self.distances2.index(max(self.distances2))
            self.positions.append(self.position2)
            
            self.area.configure(state="normal")
            message = "Calcule distance Euclidienne"+"\n"
            self.area.insert(END, message)
            self.area.configure(state="disabled")
            
            self.area.configure(state="normal")
            message = "Position dans le genome avec la max dist: "+ str(self.position2)+"\n\n"
            self.area.insert(END, message)
            self.area.configure(state="disabled")
            
            pos = 0
            for i in range(len(self.pas)):
                pos += (self.positions[i]*self.pas[i])

            self.partieGenomeDif = self.seq[pos: pos + self.fenetre2]

            self.btn4.configure(state="normal")      
            self.e4.configure(state="normal")
            self.e5.configure(state="normal")
                
            self.fig = self.myplotcode(self.subparties, self.distances2)
            self.fenetrePlot()

        except:
            mbox.showwarning("Warning", "Number must be Integer")                
    
    def onClick(self):
        if self.imp.get() == True:
            self.area.configure(state="normal")
            self.area.insert(END, self.partieGenomeDif)
            self.area.configure(state="disabled")   
        
def main():
    root = Tk()
    root.geometry("900x500+300+200")
    app = DockingProteines(root)
    root.mainloop()  


if __name__ ==  "__main__" :
    main()
