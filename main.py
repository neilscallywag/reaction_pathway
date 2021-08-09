import requests
import random
import operator


group_13 = ['B','Al', 'Ga','In']
group_14 = ['C','c','Si','Ge','Sn','Pb']
group_15 = ['N','n','P','As']
group_16 = ['O','o','S','Se']
group_17 = ['F','Cl','Br','I']



# Molecule is initialised 
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

#atom is initialised with certain properties
    #1. A.I index
    #2. Name
    #3 Valency min(2 numbers based on bonds)
    #4. Aromaticity based on if the character is lower cased or upper case
# To add:
    #1. Hybridisation to help form aromatic rings

class atom(object):
    def __init__(self,index,name,charge=0):
        self.index = index
        self.name = name
        self.valency = self.get_valency(name)
        self.charge = charge
        self.aromatic = True if name.islower() else False
    def attributes(self):
       print(self.index,self.name,self.valency,self.charge)
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
            
    
       
#idea extracted from pysmiles library 
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
    return None
print(z)
smile(z,g)
g.add_hydrogens()
print(g.graph())




            
