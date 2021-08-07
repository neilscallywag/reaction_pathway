import pprint
import requests
import random
import operator
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
pp = pprint.PrettyPrinter(indent=4)   
class molecule:
    def __init__(self):
        self.molecule = {} 
    def add_atom(self,atom):
            self.molecule[atom] = [] 
    def graph(self):
        d = {}
        for k,v in self.molecule.items():
            z = [x[0].index for x in v]          
            d[k.index] = z
        return d
             
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
    def bond_sum(self,atom):
        s = 0
        for x,y in self.molecule[atom]:
            s += y
        return s          
    def add_hydrogens(self):
        pass
    def fix_valence(self):
        pass
    def fix_charge(self):
        pass

               
class atom(object):
    def __init__(self,index,name,valency=None,charge=0):
        self.index = index
        self.name = name
        self.valency = valency
        self.charge = charge
    def attributes(self):
       print(self.index,self.name,self.valency,self.charge)
    def valency(self):
        return self.valency
    def index(self):
        return self.index
    def name(self):
        return self.name
       
#idea extracted from pysmiles library 
def tokenise(smiles):
    o = 'B C N O P S F Cl Br I * b c n o s p'.split()
    smiles = iter(z)
    queue =[]
    token = '' 
    valency_table = {"B": (3,), "C": (4,), "N": (3, 5), "O": (2,), "P": (3, 5), "S": (2, 4, 6), "F": (1,), "Cl": (1,), "Br": (1,),"I": (1,)}
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

z = random.choice(['NNccc(ccc[N+]([O-])=O)[N+]([O-])=O',
                   '[Cu+2].[O-]S(=O)(=O)[O-]',
                   'OCCc1c(C)[n+](cs1)Cc2cnc(C)nc2N',
                   'COc1cc(C=O)ccc1O',
                   'CC(=O)OC1=CC=CC=C1C(=O)O',
                   'C(C1C(C(C(C(O1)O)O)O)O)O',
                   'C(I)(I)I',
                   'CCOCC'])
g = molecule()

def smile(z,g):
    anchor = None
    index = 0
    next_bond = None
    branches= []
    default = 1
    stack = tokenise(z)
    rings = {}
    pp.pprint(stack)
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
        elif tokentype == 'be':
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
    return None
print(z)
smile(z,g)   
lol = g.graph()
lolz = nx.Graph(lol)
nx.draw_networkx(lolz, with_labels = True, node_color = "c", edge_color = "k", font_size = 8)

plt.axis('off')
plt.draw()
plt.savefig("graph.pdf")

            
