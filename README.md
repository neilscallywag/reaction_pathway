# reaction_pathway
script to solve reaction pathways for alevel questions


###Input
* Solve(initial,final,number_of_stages): 
* Take intitial and final product in the form of the chemical formula
* Take number of stages as number of steps it takes to reach the product (To limit searching through different branches of nodes)


### Backend requirements
1. Hall mark table: -
* This represents what are the prerequisities for the molecule to be considered to undergo a particular reaction. Eg
* 
| Hall Mark  | Potential reactions | Product |
| ------------- | ------------- |  ------------- | 
| If Contains Halogen(Array) && halogen bonded to carbon && Carbons_initial - Carbon_final > 0  | Dehalogenation  | C=C where C-Br bond was located |

