# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 14:28:54 2024

@author: sgsreidg
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw, Lipinski, Crippen, Descriptors, QED
from rdkit.Chem.Draw import rdMolDraw2D

suppl = Chem.SDMolSupplier('M:/Documents/3rd Year/366 project/Gold files/Frag library docking/CSVs/RunS2.sdf' , Chem.RemoveHs==False)

# create empty lists to be filled with data from the sdf file
tNames = [] # name
tMW = [] # molecular weight
tALOGP = [] # LOGP
tHBA = [] # Hydrogen Bond Acceptors
tHBD = [] # Hydrogen Bond Donors
tPSA = [] # ? 
tROTB = [] # Number of rotable bonds
tAROM = [] # number of aromatic systems
tALERTS = [] #ALERTS
tASP = [] # ASP score
mollist = [] # molecule data


for x in range (0 , len(suppl)): # create a temporary dataframe to add to the main one
    tlist = list(QED.properties(suppl[x])) #find properties of the molecule in the order: MW, ALOGP, HBA, HBD, PSA, ROTB, AROM, ALERTS
    tname = suppl[x].GetProp('_Name') # find the name
    z_pos = tname.find("Z") # find the position of the letter Z (starting position of the name)
    name = tname[z_pos:z_pos+10] # select the molecue name and not any other data stored in the name.
    
    tNames.append(name)
    tMW.append(tlist[0])
    tALOGP.append(tlist[1])
    tHBA.append(tlist[2])
    tHBD.append(tlist[3])
    tPSA.append(tlist[4])
    tROTB.append(tlist[5])
    tAROM.append(tlist[6])
    tALERTS.append(tlist[7])
    tASP.append(float(suppl[x].GetProp('Gold.ASP.Fitness')))
    # mollist.append(suppl[x])
    
mollist = [x for x in suppl]
    

"""Calculating ligand efficiency"""
LEi =[] # empty list to be filled with ligand efficiency
for mol in range (0 , len(suppl)):
    HA = rdkit.Chem.Lipinski.HeavyAtomCount(suppl[mol]) # find number of heavy atoms
    tLE = tASP[mol] / HA # calculate Approximate ligand efficiency for this molecule
    LEi.append(tLE) # add approximate ligand efficiency to list
LE = np.array(LEi) # change list into array before putting it in the dataframe

df = pd.DataFrame({'Name':np.array(tNames), # find name of molecule and create dataframe onto which the other molecules data will be added.
                    'MW': np.array(tMW),
                   'ALOGP': np.array(tALOGP),
                   'HBA': np.array(tHBA),
                   'HBD': np.array(tHBD),
                   'PSA': np.array(tPSA),
                   'ROTB': np.array(tROTB),
                   'AROM': np.array(tAROM),
                   'ALERTS': np.array(tALERTS),
                   'ASP Fitness':np.array(tASP),
                   'Ligand Efficiency':LE , 
                   'mollist': mollist})


ax1 = df.plot.scatter(x='HBA',
                      y='ASP Fitness',
                      xticks=range(0,6) ,
                      title='Effect of number of Hydrogen bond acceptors on ASP Fitness')

ax2 = df.plot.scatter(x='HBD',
                      y='ASP Fitness',
                      xticks=range(0,4) ,
                      title='Effect of number of Hydrogen bond Donors on ASP Fitness')

ax3 = df.plot.scatter(x='MW',
                      y='ASP Fitness',
                      title='Effect of number of Hydrogen bond Donors on ASP Fitness')

ax3 = df.plot.scatter(x='HBA',
                      y='Ligand Efficiency',
                      xticks=range(0,6),
                      title='HBA vs LE')

ax3 = df.plot.scatter(x='ALOGP',
                      y='ASP Fitness',
                      title='LOGP vs ASP')
