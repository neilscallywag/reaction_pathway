import requests
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


'''
    Valency Table:
    list of period table so that it can be used to determine
    valency of each atom
'''
group_13 = ['B','Al', 'Ga','In']
group_14 = ['C','c','Si','Ge','Sn','Pb']
group_15 = ['N','n','P','As']
group_16 = ['O','o','S','s','Se']
group_17 = ['F','Cl','Br','I']



'''
    Atomic massses

'''
atomic_mass = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811, 'C': 12.0107, 'N': 14.0067,
              'O': 15.9994, 'F': 18.9984032, 'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815386,
              'Si': 28.0855, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078,
              'Sc': 44.955912, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938045,
              'Fe': 55.845, 'Co': 58.933195, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.409, 'Ga': 69.723, 'Ge': 72.64,
              'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585,
              'Zr': 91.224, 'Nb': 92.90638, 'Mo': 95.94, 'Tc': 98.9063, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42,
              'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.760, 'Te': 127.6,
              'I': 126.90447, 'Xe': 131.293, 'Cs': 132.9054519, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116,
              'Pr': 140.90465, 'Nd': 144.242, 'Pm': 146.9151, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25,
              'Tb': 158.92535, 'Dy': 162.5, 'Ho': 164.93032, 'Er': 167.259, 'Tm': 168.93421, 'Yb': 173.04,
              'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.9479, 'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217,
              'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.9804,
              'Po': 208.9824, 'At': 209.9871, 'Rn': 222.0176, 'Fr': 223.0197, 'Ra': 226.0254, 'Ac': 227.0278,
              'Th': 232.03806, 'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482, 'Pu': 244.0642, 'Am': 243.0614,
              'Cm': 247.0703, 'Bk': 247.0703, 'Cf': 251.0796, 'Es': 252.0829, 'Fm': 257.0951, 'Md': 258.0951,
              'No': 259.1009, 'Lr': 262, 'Rf': 267, 'Db': 268, 'Sg': 271, 'Bh': 270, 'Hs': 269, 'Mt': 278,
              'Ds': 281, 'Rg': 281, 'Cn': 285, 'Nh': 284, 'Fl': 289, 'Mc': 289, 'Lv': 292, 'Ts': 294, 'Og': 294}

'''
    Defining an Atom:
                     An Atom for our purpose is defined as an object with the
                     following properties:
                     1. An index -> index: int
                     2. Its Chemical Symbol -> name: char
                     3. Its Valency - > valency: int
                     4. Its Charge -> Charge: undefined
                     5. If the atom is aromatic or not -> aromatic: Boolean
                     6. Its atomic mass -> ar:int

                     Extra Information:
                     Aromaticity is currently defined by lower case strings
                     according to SMILES format

'''

class atom(object):
    def __init__(self,index,name,charge=0):
        self.index = index
        self.name = name
        self.valency = self.get_valency(name)
        self.charge = charge
        self.ar = atomic_mass[name] if name in atomic_mass else 0
        self.aromatic = True if name.islower() else False
    def valency(self):
        return self.valency
    def index(self):
        return self.index
    def name(self):
        return self.name
    def aromatic(self):
        return self.aromatic
    '''
        Using the valency table above to get valency
    '''
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




    '''
        Defining a Molecule:
        Molecule is represented using a modifed adjacency list structure
        Each Node and adjacent Nodes are labelled by the Atom object


        Index Representation of Methane (CH4)(Smiles:C)

            {
          "0": {
            "1": {
              "bond-order": 1
            },
            "2": {
              "bond-order": 1
            },
            "3": {
              "bond-order": 1
            },
            "4": {
              "bond-order": 1
            }
          },
          "1": {
            "0": {
              "bond-order": 1
            }
          },
          "2": {
            "0": {
              "bond-order": 1
            }
          },
          "3": {
            "0": {
              "bond-order": 1
            }
          },
          "4": {
            "0": {
              "bond-order": 1
            }
          }
        }

        

        Object Representation of Methane:

        
        {
        <__main__.atom object at ADDR>: [(<__main__.atom object at ADDR>,1)],
        
        <__main__.atom object at ADDR>: [(<__main__.atom object at ADDR>,1)],
        
        <__main__.atom object at ADDR>: [(<__main__.atom object at ADDR>,1)],
        
        <__main__.atom object at ADDR>: [(<__main__.atom object at ADDR>,1),
                                         (<__main__.atom object at ADDR>,1),
                                         (<__main__.atom object at ADDR>,1),
                                         (<__main__.atom object at ADDR>,1)],
                                                       
        <__main__.atom object at ADDR>: [(<__main__.atom object at ADDR>,1)]
        }


    '''

class molecule:
    def __init__(self):
        self.molecule = {}
    def __iter__(self):
        for each in self.molecule.values():
            yield each

    '''
        Adding an ATOM in an empty dictionary as Key with empty list as its values.
        Indicating no adjacent atoms.

        Input is an ATOM object
        Output: NIL
    '''

    def add_atom(self,atom):
            self.molecule[atom] = []

    '''
        Printing Index Notation of the atom.
    '''
    def graph(self):
        d = {}
        for k,v in self.molecule.items():
            z = {x[0].index: {'bond-order':x[1]} for x in v}
            d[k.index] = z
        return d

    '''
        Printing the Object Representation of the Molecule
    '''
    def schema(self):
            return self.molecule

    '''
        Returns list of nodes
    '''

    def nodes(self):
        return [x.name for x in self.molecule.keys()]

    '''
        Adding a bond and the bond order between the atoms. Inputs are all in ATOM object.
        u: starting atom
        v: ending atom
        w: is the bond order. If bond order is not specified, default single bond is added
    '''
    
    def add_bond(self,u,v,w=1):
        self.molecule[u].append((v,w))
        self.molecule[v].append((u,w))
        

    '''
        Getting the number of surrounding atoms.
        It will be later used to calculate remaining valency of the attached atom so that Hydrogens can be added.

        Input is the ATOM Object
        Output is an int value.
    '''
    def get_surrounding_atoms(self,atom):
        current = len(self.molecule[atom])
        if atom in self.molecule[atom]:
            current += 1
        return current

    '''
        It is not possible for us to input an object for the methods described above.
        Hence we need to use the Index to make reference to the Atom object.

        This method returns the Atom based on its index attribute.
    '''

    def get_atom(self,index):
        k = list(self.molecule)
        key = k[index]
        return key


    '''
        Now we will start adding Hydrogens to our Molecule based on remaining valency

        We will get the valency of each ATOM currently in the Molecule.
        Next, we will calculate the sum of all the bond orders in the surrounding atoms.
        Lastly, we will minus off the sum of bond orders from the valency to get the number of hydrogens (Hydrogen Count) we can attach to that Atom.

    '''

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

    '''
        Adding the hydrogens

        1. Get the last index of the Hydrogen-Less molecule and + 1 to indicate the Atom index of the first Hydrogen
        2. Create New Atom for each Hydrogen and attach it using the methods described above

    '''
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


    '''
        Finally, Converting the SMILES formatted string to join atoms, branches, rings and generating the molecule itself.

    '''

    def smile(self,z):
        anchor = None
        index = 0
        next_bond = None
        branches= []
        default = 1
        stack = tokenise(z)
        rings = {}
        for tokentype,token in stack:
            if tokentype == 'atom':
                self.add_atom(atom(index,token))
                if anchor is not None:
                    if next_bond is None:
                        next_bond = default
                    if next_bond:
                        self.add_bond(self.get_atom(anchor),self.get_atom(index), next_bond)
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
                        self.add_bond(self.get_atom(index-1),self.get_atom(jdx),next_bond)
                        next_bond = None
                    next_bond = None
                    del rings[token]
                else:
                    rings[token] = (index-1,next_bond)     
        return self.add_hydrogens()


    '''
        Get molecular mass of the molecule
        
    '''
    def molmass(self):
        mr = 0
        for atoms in self.molecule.keys():
            mr += atoms.ar
        return mr

    '''
        Function group methods
    '''
    def exists(self):
        
        valid_fg = ['primary amine','secondary amine', 'tertiary amine',
                    'primary alcohol', 'secondary alcohol', 'tertiary alcohol',
                    'carboxylic acid', 'acyl chloride', 'ester',
                    'halogenoalkane', 'aldehyde', 'ketone']
        present = []

        '''
                Amines: Since nitrogen has a normal valence of three,
                        we can also conclude that there are two N-H bonds in primary amines
                        and one N-H bond in secondary amines.
                        In tertiary amines there are no N-H bonds.

                     
               Alcohols: Carbon needs to be bonded to an O by bond order of 1 (Single Bond).
                         The Oxygen should only be bonded to the C and another H hence len = 2
                         The other 3 peripheral groups can be either R groups or H.
                         For Primary Alcohol, C = [ H, O, H, R] & O = [ C, H ]
                         For Secondary Alcohol, C = [ H, O, R, R] & O = [ C, H ],
                         For Tertiary Alcohol, C = [ R, O, R, R] & O = [ C, H ]
         
        '''
        for atom in self.molecule.keys():

            if atom.name == 'N':
                if len(self.molecule[atom]) == 3:
                    hcount = 0
                    others = 0
                    for tuples in self.molecule[atom]:
                        if tuples[0].name == 'H':
                            hcount += 1
                        else:
                            others += 1
                    if hcount == 2 and others == 1:
                        present.append(valid_fg[0]) #primary amine
                    elif hcount == 1 and others == 2: 
                        present.append(valid_fg[1]) #secondary amine
                    elif hcount == 0 and others == 3: 
                        present.append(valid_fg[2]) #tertiary amine


            elif atom.name == 'C':
                if len(self.molecule[atom]) == 4: # C = [ R, O, R ,R]
                    hcount = 0
                    alcohol = 0
                    others = 0
                    for tuples in self.molecule[atom]:
                        if tuples[0].name == 'O' and tuples[1] == 1:
                            if len(self.molecule[tuples[0]]) == 2:# O = [C, R]
                                for tups in self.molecule[tuples[0]]:
                                    if tups[0].name == 'H':
                                        alcohol += 1
                        elif tuples[0].name == 'H':
                            hcount += 1
                        else:
                            others += 1
                    if alcohol == 1 and hcount == 2 and others == 1:   #CH2OH
                        present.append(valid_fg[3]) #primary alcohol  
                    elif alcohol == 1 and hcount == 1 and others == 2: #CHROH
                        present.append(valid_fg[4]) #secondary alcohol
                    elif alcohol == 1 and hcount == 0 and others == 3: #CR2OH
                        present.append(valid_fg[5]) #tertiary alcohol                        
            

        return present

    '''
                    End Molecule Class


                    
                    Interpreting a SMILES string. The code below is inspired by the PySMILES library.

                    Input is a SMILES string.

    '''


def tokenise(smiles):
    o = 'B C N O P S F Cl Br I * b c n o s p'.split()
    smiles = iter(smiles)
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
        elif char in '-=#$:.@':
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
            elif char == '@':
                order = 1   
            queue.append(('bond',order))
        elif char == '(':
            queue.append(('sb', '('))
        elif char == ')':
            queue.append(('eb', ')'))
        elif char.isdigit() :
            queue.append(('ring',char))
    return queue              



'''
        Fetching SMILES string from a Given IUPAC name using an API
'''

def convert(ids):
    try:
        url = f"http://cactus.nci.nih.gov/chemical/structure/{ids}/smiles"
        ans = requests.get(url).text
        return ans
    except Exception as e:
        l.fatal({traceback.print_exc()})  # haven't tested this yet lol
        




'''

Main Prompt

'''

def main():


    ids = input("Enter the IUPAC name of a molecule: ")
    print('You have entered: ', ids)
    z = convert(ids)
    print('SMILES Format: ', z)
    g = molecule()
    g.smile(z)
    a = g.graph()
    b = g.molmass()
    print('Molecular mass: ',b)
    print('functional groups present: ',g.exists())
    G = nx.Graph(a)
    pos=nx.spring_layout(G, k=1*1/np.sqrt(len(g.nodes())), iterations=999)
    nx.draw_networkx(G,pos,node_size=1, font_size=10)
    labels = nx.get_edge_attributes(G,'bond-order')
    nx.draw_networkx_edge_labels(G,pos,edge_labels=labels,font_size=5)
    plt.axis('off')
    plt.draw()
    plt.show()
	

while True:
    main()
