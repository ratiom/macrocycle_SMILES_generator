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
from collections import Counter


'''
'al': 	    {'Code': '',
                'SMILES': '',
                },
'''
#################
# Picking the amino scids for use in the library
resi_al = ["al2", "al3"] # alpha L-amino acids 
resi_ad = ["ad1"] # alpha D-amino acids
resi_bl = ["bl1"] # beta-amino acids
resi_nl = ["nl1", "nd1"] # L- and D-amino acids


resi_all = resi_al + resi_ad + resi_bl # Combine families of residues into one list

################
# Possible amino acids for use 

aminos = {
'al2': 	    {'Code': 'Gly',
                'SMILES': 'NCC(=O)O',
                },
'al3':          {'Code': 'L-Ala',
               'SMILES': 'N[C@@]([H])(C)C(=O)O'
                },
'al4':          {'Code': 'L-Arg',
                'SMILES': 'N[C@@]([H])(CCCNC(=N)N)C(=O)O'
                },
'al14': 	    {'Code': 'L-Phe',
                'SMILES': 'N[C@@]([H])(Cc1ccccc1)C(=O)O',
                },
'al21': 	    {'Code': 'L-hPhe',
                'SMILES': 'N[C@@]([H])(CCc1ccccc1)C(=O)O',
                },                
'al15': 	    {'Code': 'L-Pro',
                'SMILES': 'N1[C@@]([H])(CCC1)C(=O)O',
                },                
'ad1':          {'Code': 'D-Ala',
               'SMILES': 'N[C@]([H])(C)C(=O)O'
                },
'ad2':          {'Code': 'D-Arg',
                'SMILES': 'N[C@]([H])(CCCNC(=N)N)C(=O)O',
                },
'ad14': 	    {'Code': 'D-Phe',
                'SMILES': 'N[C@]([H])(Cc1ccccc1)C(=O)O',
                },
'bl1':          {'Code': 'b-Ala',
                 'SMILES': 'NCCC(=O)O'
                },
'nl1': 	    {'Code': 'NMe-L-Ala',
                'SMILES': 'N(C)[C@@]([H])(C)C(=O)O',
                },                
'nd1': 	    {'Code': 'NMe-D-Ala',
                'SMILES': 'N(C)[C@]([H])(C)C(=O)O',
                },
'link1':     {'Code': 'Link_SSS',
                'SMILES': 'N[C@@H](C)[C@H]5(C(=O)NC(C)(C)C)',
                },
'link2':     {'Code': 'Link_SRS',
                'SMILES': 'N[C@H](C)[C@H]5(C(=O)NC(C)(C)C)',
                },
'link3':     {'Code': 'Link_SRR',
                'SMILES': 'N[C@H](C)[C@@H]5(C(=O)NC(C)(C)C)',
                }
}

########################
#   These residues are the non-X residues in the patterns.
holding = {'P': 'al15', # L-Pro
           'F': "al14", # L-Phe
           'H': "al21", # L-homoPhe
           'f': 'ad14'  # D-Phe
           }
'''
Generate list
Give list non-X elems
Give possible lengths, i.e. 5 (PFXXX), 6 (PXFXXX)
Loop through range of lengths, to generate lists with #X = length - len(non-X)
List

li = []
ranges = [5,6]
nonX = ['F', 'f', 'H']
for i in 
'''

letters_Ff3 = ['F','f','X']
letters_Ff4 = ['F','f','X','X']
letters_Ff5 = ['F','f','X','X', 'X']
letters_FF3 = ['F','F','X']
letters_FF4 = ['F','F','X','X']
letters_FF5 = ['F','F','X','X', 'X']
letters_FH3 = ['F','H','X']
letters_FH4 = ['F','H','X','X']
letters_FH5 = ['F','H','X','X', 'X']
letters_all = [letters_Ff3, letters_Ff4, letters_Ff5, letters_FF3, letters_FF4, letters_FF5, letters_FH3, letters_FH4, letters_FH5]


def patternGen(letters_all):
    patterns = []
    for i in letters_all:
        keywords = ['P'+''.join(j) for j in itertools.product(i, repeat = len(i))]
        keywords2 = [q for q in keywords if q.count('F') == 1 and q.count('f') == 1]
        holding_count = Counter(keywords)
        #keywords2 = [q for q in keywords if q.count('F')  < 3]
        patterns = patterns + keywords2
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
            print sequence
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
        
        for v in range(1,3):
            cyclic = linear + aminos['link'+str(v)]['SMILES'] #Add on the linker SMILES
            cyclic = re.sub(r'^N1', r'N15', cyclic) #Adjust the SMILES to install the ring at Pro-N
            cycSeq = linSeq + ',' + aminos['link'  + str(v)]['Code'] #Add the linker code
            mol = (cyclic, cycSeq)
            cycles.append(mol)
    return cycles
     
patternGen(letters_all)
pepGen(resi_all, holding, aminos, patterns)
molGen(peptides)
   
sampled_cycles = random.sample(cycles, 5000) # Select random subset of the peptide list


#############
### Write list of peptides to CSV
        
os.chdir('C:\Users\Conor')
smiles_file = open('macro_smiles_random.csv', 'w')
smiles_writer = csv.writer(smiles_file)
smiles_writer.writerow(['SMILES','title'])
smiles_writer.writerows(sampled_cycles)
smiles_file.close()