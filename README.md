# reaction_pathway
script to solve reaction pathways for alevel questions


### Input
* Solve(initial,final,number_of_stages): 
* Take intitial and final product in the form of the chemical formula
* Take number of stages as number of steps it takes to reach the product (To limit searching through different branches of nodes)


### Backend requirements
1. Labelling the elements
2. Hall mark table: -
* This represents what are the prerequisities for the molecule to be considered to undergo a particular reaction. Eg
* 
| Hall Mark  | Potential reactions | Product |
| ------------- | ------------- |  ------------- | 
| If Contains Halogen(Array) && halogen bonded to carbon && Carbons_initial - Carbon_final > 0  | Dehalogenation  | C<sub>3</sub>=C<sub>4</sub> where C<sub>4</sub>-Br bond was located |






### Assumptions
* No non carbon molecule is bonded to another non carbon molecule
* No pathway where multiple intermediate products are created and we have to pick one (For now)
