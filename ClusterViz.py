import graphviz
import nltk
from collections import Counter
import re
from copy import copy


def chartrigram(string):
    liste = []
    if 3 < len(string):
        for p in range(len(string) - 2) :
            tg = string[p:p+3]
            liste.append(tg)
    return liste

def tgcharoverlap(a,b):
    tga  = chartrigram(a)
    tgb  = chartrigram(b)
    intersec = 0
    for t in tga:
        if t in tgb:
            intersec += 1
    union = len(tga) + len(tgb) - intersec
    if union == 0:
        return 0.0
    return intersec/union

def substitutionsfehler(a,b):
    error = 2
    
    if a == b:
        return 0
    elif a.lower() == b.lower():
        return  0.1
    elif len(a) < 4 and len(b) < 4:
        return 1
    else:
        return  2 - tgcharoverlap('#'+a.lower()+'#','#'+b.lower()+'#')
        

def penalty(A,B):
    return min([substitutionsfehler(a,b) for b in B for a in A])
        
        
def align(v,w,bp,n):
    empty_element = '####'
    al = []
    i,j= len(v),len(w)  
    while (i,j)  != (0,0):
        al.append((i,j)) #erst hinzufügen, sonst fehlt das letzte Token in jedem Satz
        i,j = bp[i][j]
    al = al[-1::-1]
    
    table = []
    i,j = 0,0
    for i_1,j_1 in al:
        if i_1 > i and j_1 > j:
            table.append(v[i]+w[j])
        elif i_1 > i:   
            lenb = n - len(v[i])
            table.append(v[i]+lenb*[empty_element]) #GGf. mehrere empty elemente!!
        elif j_1 > j:   
            lena = n - len(w[j])
            table.append(lena*[empty_element]+w[j]) #GGf. mehrere empty elemente!!

        i,j = i_1,j_1
    return table

def edit_distance(v, w):
    #v und w  sind jeweils listen von Varianten
    
    matrix = [[0 for j in range(len(w) + 1)] for i in range(len(v) + 1)] 
    backpointer = [[(0,0) for j in range(len(w) + 1)] for i in range(len(v) + 1)]
    for i in range(len(v)+1):
        for j in range(len(w)+1):
            if i > 0 and j > 0:
                val1 = matrix[i-1][j] + 1
                val2 = matrix[i][j-1] + 1
                val3 = matrix[i-1][j-1]  + penalty(v[i-1],w[j-1]) 
                matrix[i][j] = min(val1, val2, val3)
                if matrix[i][j] == val3:
                    backpointer[i][j] = (i-1,j-1)
                elif matrix[i][j] == val2:
                    backpointer[i][j] = (i,j-1)
                elif matrix[i][j] == val1:
                    backpointer[i][j] = (i-1,j)
            elif i > 0:
                matrix[i][j] = matrix[i-1][j] + 1
                backpointer[i][j] = (i-1,j)
            elif j > 0:
                matrix[i][j] = matrix[i][j-1] + 1
                backpointer[i][j] = (i,j-1)
            else:
                matrix[i][j] = 0 #Die erste Zelle

    ed = matrix[len(v)][len(w)]
    #alignment = align(v,w,backpointer,empty_element)
    return ed, backpointer# alignment

def aligned_sequence(sentlist):
    
    table = [[t] for t in sentlist[0]]

    n = 1
    if len(sentlist) > 1: 
        for v in sentlist[1:]: #erstes ist oben, 
            n+=1
            v1 = [[t] for t in v]
            #_, seq = edit_distance(seq, v1,'####')
            _, bp = edit_distance(table, v1)     
            table = align(table,v1,bp,n) #number of columns in alignment table
                   
    return table


class Graph:
    nodes = []
    arcs = {}
    markednodes = []
    last_id = 0

    def __init__(self):
        self.nodes = []
        self.arcs = {}
        self.last_id = 0
        
    def addnode(self,txt):
        self.last_id += 1
        self.nodes.append((self.last_id,txt))
        return self.last_id 
        
    def addarc(self,n1,n2):
        #self.arcs.append((n1,n2))
        children = self.arcs.get(n1,[])
        if n2 not in children:
            children.append(n2)
        self.arcs[n1] = children
                
    def getnode(self,nid):
        for n,txt in self.nodes:
            if n == nid:
                return txt
        return ""
        
    def children(nid):
        return self.arcs[nid]

        
    def simplify(self):
        indegree = Counter()
        for ingoing in self.arcs.values():
            indegree.update(ingoing)
        
        todo = copy(self.nodes)
        
        #for nid,txt in self.nodes:
        while len(todo) > 0:
            nid,txt = todo.pop(0)
            if nid > 1 and len(self.arcs.get(nid,[])) == 1:
                cid = self.arcs[nid][0]
                childtext = self.getnode(cid)
                if indegree[cid] == 1 and len(self.arcs.get(cid,[])) > 0: #outdegree should not be 0 to keep a seperate end-node
                    self.nodes.remove((cid,childtext))
                    if (cid,childtext) in todo:
                        todo.remove((cid,childtext))
                    if childtext == '####':
                        newtext =  txt
                    else:
                        newtext = txt + ' ' + childtext
                    self.nodes.append((nid,newtext))
                    todo.append((nid,newtext))
                    self.nodes.remove((nid,txt))  
                    clist = self.arcs[nid]
                    clist.remove(cid)
                    if cid in self.arcs:
                        clist.extend(self.arcs[cid])
                        self.arcs.pop(cid)
                    self.arcs[nid] = clist
                    
    def normalize(self,text):
        tokens = text.split()
        linelen = 0
        lbtext = ""
        for t in tokens:
            lbtext += t
            linelen += (1 + len(t))
            if linelen > 50:
                linelen = 0
                lbtext += "\n" 
            else:
                lbtext += ' ' 

        return lbtext 
            
    
    def display_graph(self,wdir,fname):
        g = graphviz.Digraph(fname, directory = wdir, strict=True) #strict damit es keine doppelten Kanten gibt 
        g.graph_attr['rankdir'] = 'TB'
        
        for nid,ntext in self.nodes:
            nid = "n" + str(nid)
            if ntext == '<START>':
                 g.node(nid,'<START>', shape='Mdiamond')
            elif  ntext == '<END>':
                g.node(nid,'<END>', shape='Msquare')
            elif int(nid[1:]) in self.markednodes:
                g.node(nid,self.normalize(ntext),fontname="times-bold")
            else:
                g.node(nid,self.normalize(ntext))
        
        for s in self.arcs:
            for t in self.arcs[s]:
                g.edge("n" + str(s),"n" + str(t))
                    
        #g.view() 
        g.render(format='png')

def table2graph(table,markedcolumn):
    G = Graph()
    startid = G.addnode('<START>')

    nrOfColumns = len(table[0])
    prevrowids = nrOfColumns * [-1]
    nodes = {}
    for txt in set(table[0]):
        if txt != '####':
            nid = G.addnode(txt)
            nodes[txt] = nid
    
    for i in range(nrOfColumns):
        txt = table[0][i]
        if txt == '####':
            prevrowids[i] = startid
        else:
            nid = nodes[txt] 
            prevrowids[i] = nid  
            G.addarc(startid,nid)
            if i == markedcolumn:
                G.markednodes.append(nid)
        
            
    for r in range(1,len(table)):
        thisrowids =  nrOfColumns * [-1]
        nodes = {}
        for txt in set(table[r]):
            if txt != '####':
                nid = G.addnode(txt)
                nodes[txt] = nid
        
        
        for i in range(nrOfColumns):
            txt = table[r][i]
            if  txt == '####':
                thisrowids[i] = prevrowids[i]
            else:
                nid = nodes[txt] 
                thisrowids[i] = nid  
                G.addarc(prevrowids[i],nid)
                if i == markedcolumn:
                    G.markednodes.append(nid)
            
        prevrowids = thisrowids
        
    endid = G.addnode('<END>')
    for nid in prevrowids: 
        G.addarc(nid,endid)
        
    return G   

def visualize(sentencelist,fname,wdir,mainsent = 0):
    tokenized_lists = [nltk.word_tokenize(s) for s in sentencelist] 
    aligned_table = aligned_sequence(tokenized_lists)
    graph = table2graph(aligned_table,mainsent)
    graph.simplify()
    #print(graph.nodes)
    #print(graph.markednodes)
    graph.display_graph(wdir,fname)
    
    
    