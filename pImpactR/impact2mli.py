import os
import numpy as np
from data import dictClass
from copy import deepcopy as copy
import MLI_getElem as getElem

__all__ = ['impact2mli']

def impact2mli(impatBeam,impactLattice,linename='impactLattice'):
  
  lattice = []
  nBend = 1
  nNL = 1
  nDr = 1
  nQuad = 1
  nThinMulti = 1
  for elem in impactLattice:
    if 'length' in elem.keys():
      L = elem.length
    type = elem.type
    if type == 'dipole':
      if elem.file_id <= 300:
        f = [getElem.sbend(name='dipole'+str(nBend),l=L,angle=elem.bending_angle,hgap=elem.pipe_radius,fint=elem.fringe_field_integration)]
      else:
        f = [getElem.sbend(name='dipole'+str(nBend),l=L,angle=elem.bending_angle,hgap=0.5*elem.pipe_radius,fint=elem.fringe_field_integration)]
      nBend = nBend+1
    elif type == 'nonlinear_insert':
      NL_k = 4*np.sin(np.pi*elem.tune_advance)**2/L
      NL_t =-elem.strength_t
      NL_c = elem.transverse_scale_c
      f = [getElem.nlinsert(name='NL'+str(nNL),zstart=0.0,zend=L,steps=elem.n_sckick*elem.n_map,zlen=L,k=NL_k,tau=NL_t,c=NL_c)]
      nNL = nNL+1
    elif type == 'drift':
      f = [getElem.drift(name='drift'+str(nDr),l=L)]
      nDr = nDr+1
    elif type == 'quad':
      f = [getElem.quadrupole(name='quad'+str(nQuad),l=L,k1=elem.B1)]
      nQuad = nQuad+1
    elif type == 'multipole_thin':
      f = [getElem.thlm(name='thlm'+str(nThinMulti),k1l=elem.KL_quad,k2l=elem.KL_sext,k3l=elem.KL_oct)]
      if elem.KL_deca!=0:
        UserWarning('KL_deca of impact multipole_thin is not implemented for conversion to MLI')
      if elem.KL_dodeca!=0:
        UserWarning('KL_dodeca of impact multipole_thin is not implemented for conversion to MLI')
      nThinMulti = nThinMulti+1
    else:
      print('Impact elem type '+type + ' is not recognized from MLI. skipping...')
      continue
    lattice=lattice+f
  line = getElem.line(name=linename,elemList=lattice)
  
  mass = impatBeam.mass*1.0e-9
  charge = impatBeam.charge
  ekinetic = impatBeam.kinetic_energy*1.0e-9
  beam = getElem.beam(mass=mass,charge=charge,ekinetic=ekinetic)
  unit = getElem.units()
  default = getElem.globaldefaults(lfrngsbend=1,tfrngsbend=1,driftexact=1)

  return [beam,unit,default]+lattice, line

  