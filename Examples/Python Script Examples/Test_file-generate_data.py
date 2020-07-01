# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 12:07:40 2020

@author: Jane
"""

#First, we generate the CSTR data with experimental noise added to the concentrations.   
#The smooth lines are the unperturbed data, and the data points are the noisy    
#measurements we use as our data points.

import numpy as np
import matplotlib.pyplot as plt

np.random.seed(40)
HMFinit = 100
LFinit = 50
Huinit = 50
taus = np.linspace (0,20,11)
Hconc = 0.1
T = [423,448]
R = 8.31446261815324*(1/1000) #units of kJ/(mol*K)

#True params
sigma = 5
A0HLF = 11.31
A0HHu = 16.69
EaHLF = 94.72
EaHHu = 141.94

try:
    del cHMF0,cHMF1,cLF0,cLF1,cHu0,cHu1
except:
    pass

size0 = len(taus)
cHMF0 = np.linspace(10,20,size0)
cHMF1 = np.linspace(10,20,size0)
cLF0 = np.linspace(10,20,size0)
cLF1 = np.linspace(10,20,size0)
cHu0 = np.linspace(10,20,size0)
cHu1 = np.linspace(10,20,size0)
total0 = np.linspace(10,20,size0)
total1 = np.linspace(10,20,size0)
kHLF0 = (10**A0HLF)*np.exp(-EaHLF/(R*T[0]))
kHLF1 = (10**A0HLF)*np.exp(-EaHLF/(R*T[1]))
kHHu0 = (10**A0HHu)*np.exp(-EaHHu/(R*T[0]))
kHHu1 = (10**A0HHu)*np.exp(-EaHHu/(R*T[1]))
#clean data
for i in range(size0):
    kHLF0 = (10**A0HLF)*np.exp(-EaHLF/(R*T[0]))
    kHLF1 = (10**A0HLF)*np.exp(-EaHLF/(R*T[1]))
    kHHu0 = (10**A0HHu)*np.exp(-EaHHu/(R*T[0]))
    kHHu1 = (10**A0HHu)*np.exp(-EaHHu/(R*T[1]))
    cHMF0[i] = HMFinit/(1+(Hconc*taus[i]*(kHLF0+kHHu0)))
    cHMF1[i] = HMFinit/(1+(Hconc*taus[i]*(kHLF1+kHHu1)))
    cLF0[i] = LFinit+(Hconc*taus[i]*kHLF0*cHMF0[i])
    cLF1[i] = LFinit+(Hconc*taus[i]*kHLF1*cHMF1[i])
    cHu0[i] = Huinit+(Hconc*taus[i]*kHHu0*cHMF0[i])
    cHu1[i] = Huinit+(Hconc*taus[i]*kHHu1*cHMF1[i])

f, ax = plt.subplots(1)
ax.plot(taus,cHMF0, label='HMF, 150C')
ax.plot(taus,cHMF1, label='HMF, 175C')
ax.plot(taus,cLF0, label='LA+FA, 150C')
ax.plot(taus,cLF1, label='LA+FA, 175C')
ax.plot(taus,cHu0, label='Humins, 150C')
ax.plot(taus,cHu1, label='Humins, 175C')
    
    
for i in range(size0):
    cHMF0[i] = cHMF0[i]+np.random.normal(0,sigma,1)
    cHMF1[i] = cHMF1[i]+np.random.normal(0,sigma,1)
    cLF0[i] = cLF0[i]+np.random.normal(0,sigma,1)
    cLF1[i] = cLF1[i]+np.random.normal(0,sigma,1)
    cHu0[i] = cHu0[i]+np.random.normal(0,sigma,1)
    cHu1[i] = cHu1[i]+np.random.normal(0,sigma,1)
 
cHMF0[0] = HMFinit
cHMF1[0] = HMFinit
cLF0[0] = LFinit
cLF1[0] = LFinit
cHu0[0] = Huinit
cHu1[0] = Huinit
    
for i in range (len(cHMF0)): 
    cHMF0[i] = max(0,cHMF0[i])
    cHMF1[i] = max(0,cHMF1[i])
    cLF0[i] = max(0,cLF0[i])
    cLF1[i] = max(0,cLF1[i])
    cHu0[i] = max(0,cHu0[i])
    cHu1[i] = max(0,cHu1[i])


ax.scatter(taus,cHMF0)
ax.scatter(taus,cHMF1)
ax.scatter(taus,cLF0)
ax.scatter(taus,cLF1)
ax.scatter(taus,cHu0)
ax.scatter(taus,cHu1)
ax.set_xlabel('Residence Time (min)')
ax.set_ylabel('Noisy CSTR Outlet Concs')
ax.set_title('Noisy CSTR Outlet Concs vs. Residence Time')
ax.legend()
plt.show()