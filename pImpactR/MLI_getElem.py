import numpy as np
import pandas as pd
from data import dictClass

#%%============================================================================
#                                global settings                               
#============================================================================== 
class cmd(dictClass)
  """
  command container
  """
  self.cmd = 'cmd'
  self.name = 'name'
  def str(self):
    f = self.name+': '+self.cmd+', '
    for k, v in self.items():
      if k not in ['name','cmd']:
        f = f + k +'='+str(v)+', '
    return f[:-2]

class beam(cmd):
  """
  beam parameter setting
  """
  def __init__(self,name='beam',particle='electron',energy='0.150511006'):
    self.cmd  = 'beam'
    self.name = name
    self.particle = particle
    self.energy = energy
  def str(self):
    f = self.name+': beam, '
    for k, v in self.items():
      if k not in ['name','cmd']:
        f = f + k +'='+str(v)+', '
    return f[:-2]
  
  
class units(dictClass):
  """
  beam parameter setting
  """
  def __init__(self,name='units',type='static'):
    self.cmd  = 'units'
    self.name = name
    self.type = type
  def str(self):
    f = self.name+': units, '
    for k, v in self.items():
      if k not in ['name','cmd']:
        f = f + k +'='+str(v)+', '
    return f[:-2]


class globaldefaults(dictClass):
  """
  beam parameter setting
  """
  def __init__(self,name='setdefaults',lfrngsbend=0,tfrngsbend=0,lfrngquad=0,tfrngquad=0,driftexact=0):
    self.cmd  = 'globaldefaults'
    self.name = name
    self.lfrngsbend = lfrngsbend
    self.tfrngsbend = tfrngsbend
    self.lfrngquad  = lfrngquad
    self.tfrngquad  = tfrngquad
    self.driftexact = driftexact
  def str(self):
    f = self.name+': globaldefaults, '
    for k, v in self.items():
      if k not in ['name','cmd']:
        f = f + k +'='+str(v)+', '
    return f[:-2]


    
#%%============================================================================
#                                  commands                               
#==============================================================================
class autotrack(dictClass):
  """
  particle tracking
  """
  def __init__(self,name='autotrack',type='symplectic',order='5'):
    self.cmd  = 'autotrack'
    self.name = name
    if type not in ['symplectic','taylor']:
      raise ValueError(type)
    self.type = type
    self.order = order
  def str(self):
    return self.name+':autotrack, type='+self.type+str(self.order)
    
    
class particledump(dictClass):
  """
  write particle data on file
  """
  def __init__(self,name='particledump',fname='rays.out',precision=9):
    self.cmd  = 'particledump'
    self.name = name
    self.fname = fname
    self.precision = precision
  def str(self):
    return self.name+':particledump, file='+self.fname+', precision='\
         + str(self.precision) +', close=false, flush=true, nunits=0' 

    
class iden(dictClass):
  """
  identity map
  """
  def __init__(self,name='clear'):
    self.name = name
  def str(self):
    return self.name+': iden'
clear = iden()

class end(dictClass):
  """
  finish MLI
  """
  def __init__(self,name='fin'):
    self.name = name
  def str(self):
    return self.name+': end'
fin = end()


class raytrace(dictClass):
  """
  read particle data
  """
  def __init__(self,name='raytrace',file1='rays.in'):
    self.name = name
    self.file1 = file1
  def str(self):
    return self.name+': raytrace, type=readonly, file1='+self.fname
  
class ptm(dictClass):
  """
  print map
  """
  def __init__(self,name='ptm',order=3,t2=0,u2=0):
    if order > 3:
      raise ValueError('order cannot be large than 3')
    self.name = name
    self.order = order
    self.t2 = t2
    self.u2 = u2
  def str(self):
    return self.name+': ptm, matrix='+str(self.order)+', poly='+str(self.order)+', t2='+str(self.t2)+', u2='+str(self.u2)
  
  
class tasm(dictClass):
  """
  calculate (chromatic, anharmonic) tune 
  """
  def __init__(self,name='tasm',iopt=1,delta=0.0,idata=1,ipmaps=0,isend=3,iwmaps=0):
    self.name = name
    self.iopt = iopt
    self.delta = delta
    self.idata = idata
    self.ipmaps = ipmaps
    self.isend = isend
    self.iwmaps = iwmaps
  def str(self):
    f = self.name+': tasm, '
    for k, v in self.items():
      if k not in ['name','cmd']:
        f = f + k +'='+str(v)+', '
    return f[:-2]
  
    
class stm(dictClass):
  """
  store map in MLI memory in location nmap
  """
  def __init__(self,name='stm',nmap=1):
    self.name = name
    self.nmap = nmap
  def str(self):
    f = self.name+': stm, '
    for k, v in self.items():
      if k not in ['name','cmd']:
        f = f + k +'='+str(v)+', '
    return f[:-2]
  
  
class gtm(dictClass):
  """
  get map from MLI memory in location nmap
    inputs: 
    nmap = MLI memory location (int)
    iopt = 1:concatenate, 2=replace map
  """
  def __init__(self,name='gtm',nmap=1,iopt=2):
    self.name = name
    self.nmap = nmap
    self.iopt = iopt
  def str(self):
    f = self.name+': gtm, '
    for k, v in self.items():
      if k not in ['name','cmd']:
        f = f + k +'='+str(v)+', '
    return f[:-2]
  
class tmo(dictClass):
  """
  output transfer map
    inputs: 
    ifile = file number s.t. fort.ifile
  """
  def __init__(self,name='tmo',ifile=17):
    self.name = name
    self.ifile = ifile
  def str(self):
    f = self.name+': tmo, '
    for k, v in self.items():
      if k not in ['name','cmd']:
        f = f + k +'='+str(v)+', '
    return f[:-2]
#%%============================================================================
#                                   lattice                               
#==============================================================================
class line(dictClass):
  """
  build line from elements
  inputs: 
    elemList = list of elements
  """
  def __init__(self,name='lattice',elemList=[]):
    self.name = name
    self.list = elemList
  def str(self):
    f = self.name+', line=('
    k=0
    for item in self.list:
      k=k+1
      if isinstance(item, str):
        f = f+ item + ' '
      else:
        try:
          f = f+item.name + ' '
        except TypeError:
          print('following is not MLI element')
          print(item)
      if k>=5:
        k=0
        f = f+ '& \n'
    f = f+')'
    return f
  
  
class tmi(dictClass):
  """
  get map from file
  inputs: 
    iopt = 1:concatenate, 2=replace map
    nopt = 1:rewind file
    iflie = read from fort.ifile 
    nskip = 
  """
  def __init__(self,name='tmi',ifile=16,iopt=2,nopt=1,nskip=0):
    self.name = name
    self.iopt = iopt
    if not isinstance(ifile,int):
      raise ValueError('ifile must be an positive integer')
    self.ifile = ifile
    self.nopt = nopt
    self.nskip = nskip
    
  def map2file(self,M=None,G=None):
    """
    write map into file fort.ifile
    inputs:
      M = matrix map (numpy array)
      G = generating polynomial (pandas DataFrame)
    """
    with open('fort.'+str(self.ifile),'w') as f:
      if isinstance(M,np.ndarray):
        for i in range(6):
          for j in range(6):
            f.write(str(i+1)+','+str(j+1)+', '+str(M[i,j])+'\n')
      if isinstance(G,pd.DataFrame):
        for i in range(len(G)):
          f.write(str(G.iloc[i].name)+', '+str(G.iloc[i].value)+'\n')
      
  def str(self):
    f = self.name+': tmi, '
    for k, v in self.items():
      if k not in ['name','cmd']:
        f = f + k +'='+str(v)+', '
    return f[:-2]
  
  
class nlinsert(dictClass):
  """
  nonlinear insert element
  inputs: 
    iopt = 1:concatenate, 2=replace map
    nopt = 1:rewind file
    iflie = read from fort.ifile 
    nskip = 
  """
  def __init__(self,name='nlinsert',zstart=0.0,zend=1.8,steps=1000,zlen=1.8,k=1.45446332708327,tau=-0.4,c=0.01):
    self.name = name
    self.zstart = zstart
    self.zend = zend
    self.steps = steps
    self.zlen = zlen
    self.k = k
    self.c = c
    self.tau = tau
    
  def str(self):
    f = self.name+': nlinsert, '
    for k, v in self.items():
      if k not in ['name','cmd']:
        f = f + k +'='+str(v)+', '
    return f[:-2]
  
  
class monitor(dictClass):
  """
  monitor
  """
  def __init__(self,name='monitor'):
    self.name = name  
  def str(self):
    f = self.name+': monitor'
    

class marker(dictClass):
  """
  marker
  """
  def __init__(self,name='marker'):
    self.name = name  
  def str(self):
    f = self.name+': marker'
    

class drift(dictClass):
  """
  quadrupole
  inputs:
    l  = length
  """
  def __init__(self,name='drift',l=0.1):
    self.name = name
    self.l = l
  def str(self):
    f = self.name+': drift, '
    for k, v in self.items():
      if k not in ['name','cmd']:
        f = f + k +'='+str(v)+', '
    return f[:-2]
    
    
class quadrupole(dictClass):
  """
  quadrupole
  inputs:
    l  = length
    k1 = gradient
  """
  def __init__(self,name='quadrupole',l=0.1,g1=None):
    self.name = name
    self.l = l
    if g1!=None:
      self.k1 = g1
  def str(self):
    f = self.name+': quadrupole, '
    for k, v in self.items():
      if k not in ['name','cmd']:
        f = f + k +'='+str(v)+', '
    return f[:-2]
    
    
class vkicker(dictClass):
  """
  vkicker
  inputs:
    l  = length
    kick = kick strengh
  """
  def __init__(self,name='vkicker',l=0.1,kick=0):
    self.name = name
    self.l = l
    self.kick = kick
  def str(self):
    f = self.name+': vkicker, '
    for k, v in self.items():
      if k not in ['name','cmd']:
        f = f + k +'='+str(v)+', '
    return f[:-2]
  
  
class sextupole(dictClass):
  """
  sextupole
  inputs:
    l  = length
    g2 = sextupole strengh
  """
  def __init__(self,name='sextupole',l=0.1,g2=None):
    self.name = name
    self.l = l
    if g2!=None:
      self.g2 = g2
  def str(self):
    f = self.name+': sextupole, '
    for k, v in self.items():
      if k not in ['name','cmd']:
        f = f + k +'='+str(v)+', '
    return f[:-2]
    
      
class sbend(dictClass):
  """
  sector dipole
  inputs:
    l = length
    angle 
    hgap
    fint = fringe field integration parameter
  """
  def __init__(self,name='sbend',l=0.1,angle=0,hgap=0,fint=0):
    self.name = name
    self.l = l
    self.angle = angle
    self.hgap = hgap
    self.fint = fint
  def str(self):
    f = self.name+': sbend, '
    for k, v in self.items():
      if k not in ['name','cmd']:
        f = f + k +'='+str(v)+', '
    return f[:-2]
      
      
class dipedge(dictClass):
  """
  dipole edge
  inputs:
    e1 = edge angle
    h 
    hgap
    fint = fringe field integration parameter
  """
  def __init__(self,name='dipedge',e1=0,h=0,hgap=0,fint=0):
    self.name = name
    self.e1 = e1
    self.h = h
    self.hgap = hgap
    self.fint = fint
  def str(self):
    f = self.name+': dipedge, '
    for k, v in self.items():
      if k not in ['name','cmd']:
        f = f + k +'='+str(v)+', '
    return f[:-2]
  
  
class srfc(dictClass):
  """
  thin rf cavity
  """
  def __init__(self,name='srfc',volts=0,freq=0):
    self.name = name
    self.volts = volts
    self.freq = freq
  def str(self):
    f = self.name+': srfc, '
    for k, v in self.items():
      if k not in ['name','cmd']:
        f = f + k +'='+str(v)+', '
    return f[:-2]