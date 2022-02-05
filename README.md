# reaction_pathway
script to help solve reaction pathways for alevel questions


### Input
* IUPAC name for a molecule 



### How it works?
1. Converting IUPAC name to SMILES string format
2. Convert SMILES string to graph data structure with Nodes as Atoms and Edges as bonds
3. Give weights to the bonds by bond order 
4. Add hydrogens based on (valency-surrounding bonds)
5. Fix aromaticity (Not much attention is given to this as most a level questions do not require u to break aromatic rings
6. Detailed explainantion of what the code does is included within the code
7. Hall mark table: -
* This represents what are the prerequisities for the molecule to be considered to undergo a particular reaction. Eg
* 
| Hall Mark  | Potential reactions | Product |
| ------------- | ------------- |  ------------- | 
| If Contains Halogen(Array) && halogen bonded to carbon && Carbons_initial - Carbon_final > 0  | Dehalogenation  | C<sub>3</sub>=C<sub>4</sub> where C<sub>4</sub>-Br bond was located |






