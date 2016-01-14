# -*- coding: utf-8 -*-
"""
Created on Sat Jan 09 09:56:19 2016

@author: Conor
"""
import itertools
import re
import csv
import os
import random



'''
'al':         {'Code': '',
                'SMILES': '',
                },
'''
#################
# Picking the amino scids for use in the library
resi_al = ["al11", "al20"] # alpha L-amino acids 
resi_ad = ["ad11", "ad20"] # alpha D-amino acids
resi_bl = ["bl1"] # beta-amino acids
resi_nl = ["nl11", "nd11", "nl20", "nd20"] # L- and D-amino acids


resi_all = resi_al + resi_ad + resi_nl # Combine families of residues into one list

################
# Possible amino acids for use 

'''
1 A, 8-G, 11 L, 14-F, 16-S, 17-T, 19-Y, 20-V
'''

aminos = {
# L-Amino Acids
'al2':         {'Code': 'Gly',
                'SMILES': 'NCC(=O)O',
                },
'al3':          {'Code': 'L-Ala',
               'SMILES': 'N[C@@]([H])(C)C(=O)O'
                },
'al4':          {'Code': 'L-Arg',
                'SMILES': 'N[C@@]([H])(CCCNC(=N)N)C(=O)O'
                },
'al11':         {'Code': 'L-Leu',
                 'Abb': 'L',
                'SMILES': 'N[C@@]([H])(CC(C)C)C(=O)O',
                },
'al14':         {'Code': 'L-Phe',
                'SMILES': 'N[C@@]([H])(Cc1ccccc1)C(=O)O',
                },
'al16':         {'Code': 'L-Ser',
                 'Abb': 'S',
                'SMILES': 'N[C@@]([H])(CO)C(=O)O',
                },
'al17':         {'Code': 'L-Thr',
                 'Abb': 'T',
                'SMILES': 'N[C@@]([H])([C@]([H])(O)C)C(=O)O',
                },
'al18':         {'Code': 'L-Tyr',
                 'Abb': 'Y',
                'SMILES': 'N[C@@]([H])(Cc1ccc(O)cc1)C(=O)O',
                },
'al20':         {'Code': 'L-Val',
                 'Abb': 'V',
                'SMILES': 'N[C@@]([H])(C(C)C)C(=O)O',
                },
'al21':         {'Code': 'L-hPhe',
                'SMILES': 'N[C@@]([H])(CCc1ccccc1)C(=O)O',
                },                
'al15':         {'Code': 'L-Pro',
                'SMILES': 'N1[C@@]([H])(CCC1)C(=O)O',
                },                

# D-Amino Acids
'ad1':          {'Code': 'D-Ala',
               'SMILES': 'N[C@]([H])(C)C(=O)O'
                },
'ad2':          {'Code': 'D-Arg',
                'SMILES': 'N[C@]([H])(CCCNC(=N)N)C(=O)O',
                },
'ad11':         {'Code': 'D-Leu',
                 'Abb': 'l',
                'SMILES': 'N[C@]([H])(CC(C)C)C(=O)O',
                },
'ad14':         {'Code': 'D-Phe',
                'SMILES': 'N[C@]([H])(Cc1ccccc1)C(=O)O',
                },
'ad15':         {'Code': 'D-Pro',
                'SMILES': 'N1[C@]([H])(CCC1)C(=O)O',
                },                
'ad16':         {'Code': 'D-Ser',
                 'Abb': 's',
                'SMILES': 'N[C@]([H])(CO)C(=O)O',
                },
'ad17':         {'Code': 'D-Thr',
                 'Abb': 't',
                'SMILES': 'N[C@]([H])([C@@]([H])(O)C)C(=O)O',
                },
'ad18':         {'Code': 'D-Tyr',
                 'Abb': 'y',
                'SMILES': 'N[C@]([H])(Cc1ccc(O)cc1)C(=O)O',
                },
'ad20':         {'Code': 'D-Val',
                 'Abb': 'v',
                'SMILES': 'N[C@]([H])(C(C)C)C(=O)O',
                },

#Beta Amino Acids
'bl1':          {'Code': 'b-Ala',
                 'SMILES': 'NCCC(=O)O'
                },
# N-Me Amino Acids
'nl1':         {'Code': 'NMe-L-Ala',
                'SMILES': 'N(C)[C@@]([H])(C)C(=O)O',
                },                
'nl11':         {'Code': 'NMe-L-Leu',
                 'Abb': '-NMe-L-Leu-',
                'SMILES': 'N(C)[C@@]([H])(CC(C)C)C(=O)O',
                },
'nl20':         {'Code': 'NMe-L-Val',
                 'Abb': '-NMe-L-Val-',
                'SMILES': 'N(C)[C@@]([H])(C(C)C)C(=O)O',
                },
'nd1':         {'Code': 'NMe-D-Ala',
                'SMILES': 'N(C)[C@]([H])(C)C(=O)O',
                },
'nd11':         {'Code': 'NMe-D-Leu',
                 'Abb': '-NMe-D-Leu-',
                'SMILES': 'N(C)[C@]([H])(CC(C)C)C(=O)O',
                },
'nd20':         {'Code': 'NMe-L-Val',
                 'Abb': '-NMe-L-Val-',
                'SMILES': 'N(C)[C@]([H])(C(C)C)C(=O)O',
                },

#Linkers
'link1':        {'Code': 'Link_SSS',
                'SMILES': 'N[C@@H](C)[C@H]5(C(=O)NC(C)(C)C)',
                },
'link2':        {'Code': 'Link_SRS',
                'SMILES': 'N[C@H](C)[C@H]5(C(=O)NC(C)(C)C)',
                },
'link3':        {'Code': 'Link_SRR',
                'SMILES': 'N[C@H](C)[C@@H]5(C(=O)NC(C)(C)C)',
                },
'link4':        {'Code': 'Link_SSR',
                'SMILES': 'N[C@@H](C)[C@@H]5(C(=O)NC(C)(C)C)',
                }
}

########################
#   These residues are the non-X residues in the patterns.
holding = {'P': 'al15', # L-Pro
           'p': "ad15", # D-Phe
           'S': "al16", # L-Ser
           's': 'ad16', # D-Ser
           'T': 'al17', # L-Thr
           't': 'D-Thr', # D-Thr
           'Y': 'L-Tyr', # L-Tyr
           'y': 'D-Tyr', # D-Tyr
           }

lengths = [5] # What the desired sequences should look like
                # i.e. PFXX is 4; PFXXX is 5 etc.

nonX = ['S', 's', 'T', 't', 'Y', 'y']  # What static residues will be appearing

how_many_NonX = 1
                        
def patternGen2(a, b):
    letters_all = []
    for residues in itertools.product(nonX, repeat = how_many_NonX):
        for i in lengths:
            pep=list(residues)
            li = []
            li=['X'] * (i-1 -(how_many_NonX))
            pep.extend(li)
            letters_all.append(pep)
            #print pep
            del pep
    return letters_all

    
##########
## To do: Need to construct some mechanism whereby certain positions can only hold certain residue
##########    
def patternGen(a):
    patterns = []
    for i in letters_all:
        keywords = ['P'+''.join(j) for j in itertools.product(i, repeat = len(i))]
        
        for y in keywords:
            
            check_nonX = len([t for t in y if t in nonX])
            
            if check_nonX != how_many_NonX:
                pass
            elif y in patterns:
                pass
            else:
               patterns.append(y)
             
    return patterns


def pepGen(resi_all, holding, aminos, patterns):
    peptides = []
    for k in patterns:
        for pep in itertools.product(resi_all, repeat=k.count('X')):
            pep = list(pep)
            outpep = []
            for resi in k:
                if  resi != 'X':
                    outpep.append(holding[resi])
                else:
                    outpep.append(pep.pop(0))
            
            sequence = [aminos[w]['Code'] for w in outpep]
            sequence = ','.join(sequence)
            #print sequence
            if len([i for i in outpep if i in resi_nl]) < 3 and len([i for i in outpep if i in resi_bl]) < 3:
                peptides.append((outpep,sequence))
            else:
                pass
    return peptides
    
    

def molGen(peptides):
    cycles = []
    for peptide,linSeq in peptides:
        linear_smiles = []
        for j in peptide:
            linear_smiles.append(aminos[j]['SMILES'][:-1])
        linear = ''.join(linear_smiles)
        
        for v in range(1,5):
            cyclic = linear + aminos['link'+str(v)]['SMILES'] #Add on the linker SMILES
            cyclic = re.sub(r'^N1', r'N15', cyclic) #Adjust the SMILES to install the ring at Pro-N
            cycSeq = linSeq + ',' + aminos['link'  + str(v)]['Code'] #Add the linker code
            mol = (cyclic, cycSeq)
            cycles.append(mol)
    return cycles


letters_all = patternGen2(lengths, nonX)


     
patterns = patternGen(letters_all)

peptides = pepGen(resi_all, holding, aminos, patterns)

cycles = molGen(peptides)
   
sampled_cycles = random.sample(cycles, len(cycles)/10 if len(cycles)/10 <= 5000 else 5000) # Select random subset of the peptide list


#############
### Write list of peptides to CSV
        
os.chdir('C:\Users\Conor')
smiles_file = open('macro_smiles_random.csv', 'w')
smiles_writer = csv.writer(smiles_file)
smiles_writer.writerow(['SMILES','title'])
smiles_writer.writerows(sampled_cycles)
smiles_file.close()