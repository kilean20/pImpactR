# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 14:29:47 2016

@author: kilean
"""

import pIMPACT as pm

#%%
beam = pm.getBeam()
beam['nCore_x']=1
beam['nCore_y']=1
beam['energy'] = 300.0E6     # MeV
beam['mass'] = pm.util.Mp   # proton mass MeV/c^2
beam['current'] = 0.02       # current averaged over bunch length ( bunch length is 1/frequency ). Set to 0 if no spage-charge
beam['n_particles'] = 4000   # number of macro particles in beam
beam['standard output']=2
beam['mesh_x']=32
beam['mesh_y']=32
beam['mesh_z']=32
#%%
drift0 = pm.getElem('drift') # get default drift
drift0['length'] = 0.2 # meter
drift0['n_sckick']= 2  # set number of space-charge kicks

drift1 = pm.getElem('drift')
drift1['length'] = 0.4
drift1['n_sckick'] = 4  
#%%
qf = pm.getElem('quad')
qd = pm.getElem('quad')
#%%
FODO = [drift0, qf, drift1, qd, drift0]

FODO[1]['B1'] =  15.0  # python list argument start from 0
FODO[3]['B1'] = -15.0

for i in range(len(FODO)):
    if FODO[i]['type']=='quad':
        FODO[i]['n_sckick']=5 # set number of space-charge kicks for quads
#%%
lattice = FODO + FODO + FODO
#%%
pm.writeIMPACT('test.in',beam,lattice)
#%%
betx, bety, betz = 8.0, 8.0, 25.0 # beta-function  [meter,meter,deg/MeV]  see impact.py
alfx, alfy, alfz = -1.0, 1.0, 0.0 # alpha-function [rad, rad, rad]
enx, eny, enz = 0.25E-6, 0.25E-6, 0.25 # normalized emittance  [mm-mrad, mm-mrad, deg-MeV]
pm.twiss2beam(beam,betx,alfx,enx,bety,alfy,eny,betz,alfz,enz) # use twiss parameters to define initial beam
pm.writeIMPACT('test.in',beam,lattice)

#%%
pm.run()#(nCore = 16)
#%%
# read reference orbit 
#pm.readReferenceOrbit()  # z[meter] , phase advance [degree], gamma, kinetic energy [MeV], beta 
#%%
print pm.readOpticsAt(0,'x') # read (betx, alfx, enx ) at the begining of the lattice
print pm.readOpticsAt(-1,'x') # read (betx, alfx, enx ) at end of the lattice
print pm.readOpticsAt(2,'x') # read (betx, alfx, enx ) at the i(=2)-th location 

Z_at_2 = pm.readReferenceOrbitAt(2)[0]
print 'z at i=2 is', Z_at_2   # i(=2)-th location [meter]

#%%
#pm.readBeamSize('x') # read rms beam envelope [meter]
pm.plot.rms(3,halo=95)

#%%
pm.plot.maxAmplitude()


Optics = pm.readOptics('y')
import matplotlib.pyplot as plt
plt.figure()
plt.plot(Optics[:,-1],Optics[:,0])
plt.savefig('beta.png')
plt.figure()
plt.plot(Optics[:,-1],Optics[:,-2])
plt.savefig('phadvance.png')
