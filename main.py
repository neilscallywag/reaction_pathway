import re as re

def synthesise(initial,final,steps):
    s = sanitize(initial)
    e = sanitize(final)
    r = nu_sub(s,e,steps) 
    return print(r)



def valency(i):
    pass

def functional_groups(i):
    pass
"""
def sanitize(a):
    node = re.findall("C",a)
    attachments = re.split("C", a)
    attachments =list(filter(None,attachments))
    molecule = [dict(zip(node,(a,b))) for a,b in zip(node,attachments)]
    for i in molecule:
        j = "".join(i.values())
        k = "".join(i.keys())
        l = re.sub("\d","",j)
        n = re.findall("\d",j)
        t = re.sub(j,l,j)
        r = ["Br","Cl","I","F","At","Ts","OH"]
        w =[]
        w = [z for z in t if z != ""]*int(n[0])
        for hx in r:
            stack = []
            if hx in t:
                t=t.split(hx)
                w = [z for z in t if z != ""]*int(n[0])
                stack.append(hx)
                w.append(stack.pop())   
        u = {k:list(w)}
        i.update(u)
    return molecule

"""

def text_num_split(item): #Written by: mime-r
    for index, letter in enumerate(item, 0):
        if letter.isdigit():
            return [item[:index],item[index:]]
        
def sanitize(formula): # Written by: mime-r
    # for example if formula = "CH3CH2Br"
    split_1 = re.findall('[A-Z][^A-Z]*', formula)
    # ['C', 'H3', 'C', 'H2', 'Br']
    
    # remove Cs
    split_2 = []
    count_C = split_1.count("C")
    for count in range(count_C):
        split_2.append({})
        
    n = -1
    for element in split_1:
        if element == "C":
            n += 1
            split_2[n] = ['C', []]
        else:
            r = text_num_split(element)
            if r == None:
                split_2[n][1].append(element)
            else:
                for i in range(int(r[1])):
                    split_2[n][1].append(r[0])
    # [['C', ['H', 'H', 'H']], ['C', ['H', 'H', 'Br']]]
    # now convert it into dict
    dict_1 = []
    for element in split_2:
        dict_1.append({element[0]: element[1]})
    return dict_1

def nu_sub(i,f,steps_left):
    #need to modify the input product to match intermediate
    #need to add in stack reagent and conditions
    halogens = ["Br","Cl","I","F","At","Ts"]
    for x in halogens:
        for j in i:
            n = "".join(str(y) for y in j.values())
            z = re.findall(x,n)
            if(z):
                for k in f:
                    b = "".join(str(y) for y in k.values())
                    if 'OH' in b:
                        if len(i) == len(f):
                            steps_left -= 1
                            return True , steps_left
    return False

           
#synthesise('CH3C(=O)I','CH3CH2OH',1) #error Molecule cannot undergo nu sub
#synthesise('CH3CI','CH3CH2OH',1) #error - carbon valency
synthesise('CH3CH2Br','CH3CH2OH',1) #error - Recognises the C in Cl as a carbon node


