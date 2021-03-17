import os
import numpy as np
from data import dictClass
from data import beam as _beam
from util import Me
from util import Mp
from util import cLight
from copy import deepcopy as copy
from impactIO import getElem


__all__ = ['beam','latice']
_beam=_beam()
def beam(elemList=[]):
  #----beam parameters----
  f = _beam
  f.distribution.distribution_type = 'ReadFile'
  for item in elemList:
    if item.elem == 'beam':
      # find mass
      if 'mass' in item:
        f.mass = item.mass*1.0e9
      elif 'particle' in item:
        if item.particle == 'electron':
          f.mass = Me
        elif item.particle == 'proton':
          f.mass = Mp
        else:
          print('unsupported particle type: '+item.partice+'. Manual setting is required for mass,kinetic_energy and charge.')
          f.mass = Mp
      else:
        print('Cannot determine particle mass. Manual setting is required for mass,kinetic_energy and charge.')
        f.mass = Mp
      # find energy
      if 'energy' in item:
        f.kinetic_energy = item.energy*1.0e9 - f.mass
      elif 'ekinetic' in item:
        f.kinetic_energy = item.ekinetic*1.0e9
      else:
        print('Cannot determine reference kinetic energy. Manual setting for kinetic_energy is required.')
        f.kinetic_energy = 1.0e9
      # find charge
      if 'charge' in item:
        f.charge = item.charge
      elif 'particle' in item:
        if item.particle == 'electron':
          f.charge = -1.0
        elif item.particle == 'proton':
          f.charge = 1.0
        else:
          print('unsupported particle type: '+item.partice+'. Manual setting is required for particle charge.')
          f.charge = 1.0
      else:
        print('Cannot determine particle charge. Manual setting for particle charge.')
        f.charge = 1.0
      f.multi_charge.q_m[0] = f.charge/f.mass
      break
  return f

def lattice(elemList=[],MLIline=[],impactBeam=None):
  #----build lattice----
  brho = None
  if impactBeam != None:
    g0 = impactBeam.kinetic_energy/impactBeam.mass + 1.0
    b0 = np.sqrt((g0+1.0)*(g0-1.0))/g0
    brho=g0*b0*impactBeam.mass/cLight
  fList = []
  for name in MLIline.list:
    f=None
    for item in elemList:
      if item.name==name:
        elem=item.elem
        if   elem == 'sbend':
          f = getElem('dipole')
          f.file_id = 350
          f.length = item.l
          f.pipe_radius = 2.0*item.hgap
          if f.pipe_radius == 0:
            f.pipe_radius = 1.0
          f.fringe_field_integration = item.fint
          f.bending_angle = item.angle
          f.entrance_angle = item.e1
          f.exit_angle = item.e2
          f=[f]
        elif elem == 'nlinsert':
          f = getElem('nonlinear_insert')
          f.length = item.zlen
          f.strength_t = -item.tau
          f.transverse_scale_c = item.c
          f.tune_advance = np.arcsin(np.sqrt(0.25*item.k*f.length))/np.pi
          f=[f]
        elif elem == 'drift':
          f = getElem('drift')
          f.length = item.l
          f=[f]
        elif elem == 'quadrupole':
          f = getElem('quad')
          f.length = item.l
          if 'k1' in item.keys():
            f.Kx = item.k1
          elif 'g1' in item.keys():
            if brho==None:
              raise ValueError('impactBeam is required to convert quadrupole unit from MLI g1 to impact B1')
            else:
              f.Kx = item.g1*brho
          else:
            print('quadrupole strength not found in the following MLI element:')
            print(item)
            f.Kx = 0.0
          f=[f]
        elif elem == 'thlm':
          f = getElem('multipole_thin')
          if 'k1l' in item.keys():
            f.KL_quad= item.k1l
          if 'k2l' in item.keys():
            f.KL_sext= item.k2l
          if 'k3l' in item.keys():
            f.KL_oct = item.k3l
          f=[f]
        elif elem == 'sextupole':
          L = item.l
          N = int(np.ceil(L/0.02))
          print('MLI sextupole of length '+str(L)+' found. Converting to '+str(N)+' drift-kick-drift thin element')
          kick = getElem('multipole_thin')
          drif = getElem('drift')
          if 'k2' in item.keys():
            kick.KL_sext = item.k2*L/N
          elif 'g2' in item.keys():
            if brho==None:
              raise ValueError('impactBeam is required to convert sextupole unit from MLI g2 to impact k2l')
            else:
              kick.KL_sext = item.g2*brho*L/N
          else:
            print('sextupole strength not found in the following MLI element:')
            print(item)
          drif.length = 0.5*L/N
          f=[drif,kick,drif]*N
        elif elem == 'vkicker':
          print(elem + ' is not recognized. replacing with drift...')
          f = getElem('drift')
          f.length = item.l
          f=[f]
        elif elem == 'hkicker':
          print(elem + ' is not recognized. replacing with drift...')
          f = getElem('drift')
          f.length = item.l
          f=[f]
        else:
          print(elem + ' is not recognized. skipping...')
        break
    if f!= None:
      fList=fList+f
    for elem in fList:
      if 'length' in elem:
        elem.n_sckick = int(np.ceil(elem.length*25))
  return fList

  