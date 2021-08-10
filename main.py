import requests
import random
import operator
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json
import yaml
from collections import defaultdict
from pprint import pprint

group_13 = ['B','Al', 'Ga','In']
group_14 = ['C','c','Si','Ge','Sn','Pb']
group_15 = ['N','n','P','As']
group_16 = ['O','o','S','s','Se']
group_17 = ['F','Cl','Br','I']
class molecule:
    def __init__(self):
        self.molecule = {}
    def __iter__(self):
        for each in self.molecule.values():
            yield each    
    def add_atom(self,atom):
            self.molecule[atom] = [] 
    def graph(self):
        d = {}
        for k,v in self.molecule.items():
            z = {x[0].index: {'bond-order':x[1]} for x in v}
            d[k.index] = z
        return d
    def schema(self):
        return self.molecule         
    def add_bond(self,u,v,w=1):
        self.molecule[u].append((v,w))
        self.molecule[v].append((u,w))
    def nodes(self):
        return [x.attributes() for x in self.molecule.keys()]
    def bonds(self):
        for keys in self.molecule:
            pp.pprint((keys, self.molecule[keys]))
    def get_surrounding_atoms(self,atom):
        current = len(self.molecule[atom])
        if atom in self.molecule[atom]:
            current += 1
        return current
    def get_atom(self,index):
        k = list(self.molecule)
        key = k[index]
        return key
    def get_valence(self,index):
        k = list(self.molecule)
        key = k[index]
        return key.valency        
    def bond_sum(self,atom):
        s = 0
        for x,y in self.molecule[atom]:
            s += y
        return s
    def HC(self,atom):
        s = self.bond_sum(atom)
        return atom.valency - s
    def add_hydrogens(self):
        last_atom_index = list(self.molecule)[-1]
        counter = last_atom_index.index + 1
        f = self.molecule.copy()
        for k in f.keys():
            hydrogen_count = self.HC(k)
            for hydrogen_index in range(hydrogen_count):
                self.add_atom(atom(counter,'H'))
                self.add_bond(self.get_atom(k.index),self.get_atom(counter), 1)
                counter +=1        
    def fix_valence(self):
        pass
    def fix_charge(self): #currently charges are not included 
        pass
    def fix_aromatic_bonds(self): #currently all aromatic bonds are either order 1 or 2
        pass
    def carboxyl(self) -> bool:
        c = 0
        for atom in self.molecule.keys():
            for value in self.molecule[atom]:
                print(value[0].name , value[1])
                if (value[0].name == 'O' and value[1] == 2 and (atom.name == 'C' or atom.name == 'c')) or (value[0].name == 'C' and value[1] == 2 and (atom.name == 'O')):
                    c +=1
        if c == 2:
            return print('y')
        else:
            return print('n')
                    
                    
                

        #for atom in self.molecule:
         #   if atom.name == "C":
          #      for values in atom:
           #         if values[0].name == "O" and values[1] == 2:
            #            return True
             #       else:
              #          return False
                        


class atom(object):
    def __init__(self,index,name,charge=0):
        self.index = index
        self.name = name
        self.valency = self.get_valency(name)
        self.charge = charge
        self.aromatic = True if name.islower() else False
    def attributes(self):
       return None
    def valency(self):
        return self.valency
    def index(self):
        return self.index
    def name(self):
        return self.name
    def aromatic(self):
        return self.aromatic
    def get_valency(self,name):
        if name in group_13:
            return 3
        elif name in group_14:
            return 4
        elif name in group_15:
            return 3
        elif name in group_16:
            return 2
        elif name in group_17:
            return 1
        else:
            return 0
def tokenise(smiles):
    o = 'B C N O P S F Cl Br I * b c n o s p'.split()
    smiles = iter(z)
    queue =[]
    token = '' 
    peek = None
    while True:
        char = peek if peek else next(smiles,'')
        peek = None
        if not char:
            break
        if char == '[':
            token = char
            for char in smiles:
                token += char
                if char == ']':
                    break
            queue.append(('atom',token))
        elif char in o:
            peek = next(smiles,'')
            if char + peek in o:
                queue.append(('atom',char + peek))
                peek = None
            else:
                queue.append(('atom',char))
        elif char in '-=#$:.':
            if char == '-':
                order = 1
            elif char == '=':
                order = 2
            elif char == '#':
                order = 3
            elif char == '$':
                order = 4
            elif char == ':':
                order = 1.5
            elif char == '.':
                order = 0
            queue.append(('bond',order))
        elif char == '(':
            queue.append(('sb', '('))
        elif char == ')':
            queue.append(('eb', ')'))
        elif char.isdigit() :
            queue.append(('ring',char))
    return queue                                      
def convert(ids):
    try:
        url = f"http://cactus.nci.nih.gov/chemical/structure/{ids}/smiles"
        ans = requests.get(url).text
        return ans
    except Exception as e:
        l.fatal({traceback.print_exc()})  # haven't tested this yet lol
def smile(z,g):
    anchor = None
    index = 0
    next_bond = None
    branches= []
    default = 1
    stack = tokenise(z)
    rings = {}
    for tokentype,token in stack:
        if tokentype == 'atom':
            g.add_atom(atom(index,token))
            if anchor is not None:
                if next_bond is None:
                    next_bond = default
                if next_bond:
                    g.add_bond(g.get_atom(anchor),g.get_atom(index), next_bond)
                next_bond = None
            anchor = index
            index += 1
        elif tokentype == 'sb':
            branches.append(anchor)
        elif tokentype == 'eb':
            anchor = branches.pop()
        elif tokentype == 'bond':
            next_bond = token
        elif tokentype == 'ring':
            if token in rings:
                jdx, order = rings[token]
                if next_bond is None and order is None:  
                    next_bond = 1
                elif order is None:
                    next_bond = 1
                elif next_bond is None:
                    next_bond = order
                elif next_bond != order:
                    print("error")
                    break
                elif int(index - 1) == jdx:
                    print("error2")
                if next_bond is not None:
                    g.add_bond(g.get_atom(index-1),g.get_atom(jdx),next_bond)
                    next_bond = None
                next_bond = None
                del rings[token]
            else:
                rings[token] = (index-1,next_bond)     
    return g.add_hydrogens()





z = 'CC(=O)OC'

'''
random.choice(['NNccc(ccc[N+]([O-])=O)[N+]([O-])=O',
                   '[Cu+2].[O-]S(=O)(=O)[O-]',
                   'OCCc1c(C)[n+](cs1)Cc2cnc(C)nc2N',
                   'COc1cc(C=O)ccc1O',
                   'CC(=O)OC1=CC=CC=C1C(=O)O',
                   'C(C1C(C(C(C(O1)O)O)O)O)O',
                   'C(I)(I)I',
                   'CCOCC',
                   'c1ccccc1', # need to kekulize this benzene. Current output results in all the bonds being order 1
                   'C1CCCCC1',
                   'C#N',
                   'NC(Cl)(Br)C(=O)O', #canonical
                   'O=C(O)C(N)(Br)Cl' #non canonical
                   ])

                '''
g = molecule()
print(z)
smile(z,g)
lol = g.schema()
pprint(lol)
g.carboxyl()

'''
lol = g.graph()
zzz = g.schema()
print(json.dumps(lol, indent=2, default=str))
pprint(zzz)
G = nx.Graph(lol)
pos=nx.spring_layout(G, k=0.3*1/np.sqrt(len(g.nodes())), iterations=500) 
nx.draw_networkx(G,pos)
labels = nx.get_edge_attributes(G,'bond-order')
nx.draw_networkx_edge_labels(G,pos,edge_labels=labels)
plt.axis('off')
plt.draw()
plt.savefig("graph.pdf")
'''
