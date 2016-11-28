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

if version == 3:

	#### pour python 3
    import tkinter as tkinter
    from tkinter import Tk, Text, BOTH, W, N, E, S, END, CENTER
    from tkinter.ttk import Frame, Button, Label, Style
    import tkinter.filedialog as filedialog
    from tkinter.scrolledtext import ScrolledText  
    from tkinter.font import Font
    from tkinter.messagebox import Message
    import tkinter.messagebox as mbox



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

        Label(self, text="Get Seq").grid(sticky=E, row=0, padx =5, pady=5)
        Label(self, text="Kmers").grid(sticky=E, row=1, padx =5, pady=5)

        self.btn = Button(self, text="Open File", command = self.onOpen)
        self.btn.grid(sticky=W, row=0, column=1, padx=5, pady=5)
        self.e1 = tkinter.Entry(self, state="disabled")
        self.e1.configure(width=5)
        self.e1.grid(sticky=W, row=1, column=1, padx =5, pady=5)
        
        self.btn2 = Button(self, text="Calcule", state="disabled", command = self.onCalcule)
        self.btn2.grid(sticky=W, row=2, column=1, padx=5, pady=5)
        
        lbl = Label(self, text="Sequence")
        lbl.grid(sticky=W, row=3, column=0, padx=5)
        
        lbl2 = Label(self, text="List Kmers")
        lbl2.grid(sticky=W, row=0, column=2, padx=5)
        
        self.area = ScrolledText(self, wrap=tkinter.WORD)
        self.area.configure(state="disabled", font=self.customFont, height= 20, width=50)
        self.area.grid(padx =5, pady=5, row=4, columnspan = 2)
        
        self.listbox = tkinter.Listbox(self)
        self.listbox.config(height= 25, width=20)
        self.listbox.grid(pady=5, row=1, column = 2, rowspan=5)

    def onOpen(self):
        ftypes = [('Fasta files', '*.fasta')]
        dlg = filedialog.Open (self, filetypes = ftypes)
        fl = dlg.show()
        self.seq = lireSeq(fl)
        self.area.configure(state="normal")        
        self.area.insert(END, seq)
        self.area.configure(state="disabled") 
        self.btn2.configure(state="normal")      
        self.e1.configure(state="normal")
        self.btn.configure(state="disabled") 

    def onCalcule(self):
        k = self.e1.get()
        try:
            self.kmer = int(k)
            self.e1.configure(state="disabled")
            self.btn2.configure(state="disabled")
            self.dic = find_kmers(self.seq, self.kmer)
            cont = 0
            for key in self.dic:
                string = '{} {:>10}'.format( key, str(self.dic[key]))
                print (string)
                self.listbox.insert(cont, string)
                cont += 1

        except:
            mbox.showwarning("Warning", "Number must be Integer")                
        
def main():
    root = Tk()
    root.geometry("600x500+300+200")
    app = DockingProteines(root)
    root.mainloop()  


if __name__ ==  "__main__" :
    main()  
