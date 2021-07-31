import re as re

def synthesise(initial,final,steps):
    s = sanitize(initial)
    e = sanitize(final)
    r = nu_sub(s,e,steps) 
    return print(r)
    
def sanitize(a):
    node = re.findall("C",a)
    attachments = re.split("C", a)
    attachments =list(filter(None,attachments))
    molecule = [dict(zip(node,(a,b))) for a,b in zip(node,attachments)]
    return molecule

def nu_sub(i,f,steps_left):
    #need to modify the input product to match intermediate
    #need to add in stack reagent and conditions
    halogens = ["Br","Cl","I","F","At","Ts"]
    for x in halogens:
        for j in i:
            n = "".join(j.values())
            z = re.findall(x,n)
            if(z):
                for k in f:
                    b = "".join(k.values())
                    if 'H2OH' in b:
                        if len(i) == len(f):
                            steps_left -= 1
                            return True , steps_left
    return False
                
synthesise('CH3C(=O)I','CH3CH2OH',1) #error Molecule cannot undergo nu sub
synthesise('CH3CI','CH3CH2OH',1) #error - carbon valency
synthesise('CH3CH2Cl','CH3CH2OH',1) #error - Recognises the C in Cl as a carbon node
