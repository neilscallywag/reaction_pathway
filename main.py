import requests
import random
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import json


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
    Defining an Atom:
                     An Atom for our purpose is defined as an object with the
                     following properties:
                     1. An index -> index: int
                     2. Its Chemical Symbol -> name: char
                     3. Its Valency - > valency: int
                     4. Its Charge -> Charge: undefined
                     5. If the atom is aromatic or not -> aromatic: Boolean

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


    ids = input("Enter the IUPAC name of a molecule")
    print('You have entered: ', ids)
    z = convert(ids)
    print('SMILES Format: ', z)
    g = molecule()
    g.smile(z)
    a = g.graph() #index notation
    print('Your molecule is saved as graph.pdf')
    print(json.dumps(a, indent=2, default=str)) #printing the graph notation in JSON format
    G = nx.Graph(a)
    pos=nx.spring_layout(G, k=0.3*1/np.sqrt(len(g.nodes())), iterations=500)
    nx.draw_networkx(G,pos)
    labels = nx.get_edge_attributes(G,'bond-order')
    nx.draw_networkx_edge_labels(G,pos,edge_labels=labels)
    plt.axis('off')
    plt.draw()
    plt.savefig("graph.pdf")
	

if __name__ == "__main__":
    main()
