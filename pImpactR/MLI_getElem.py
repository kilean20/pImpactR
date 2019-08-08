import numpy as np
import pandas as pd


#%%============================================================================
#                               data class                           
#============================================================================== 
def defaultKeyVal(d,k,v):
  if k in d.keys():
    return d[k]
  else:
    return v
    

class dictClass(dict):
  """ 
  This class is essentially a subclass of dict
  with attribute accessors, one can see which attributes are available
  using the `keys()` method.
  """
  def __dir__(self):
      return self.keys()
    
  def __getattr__(self, name):
    try:
      return self[name]
    except KeyError:
      raise AttributeError(name)
  if dict==None:
    __setattr__ = {}.__setitem__
    __delattr__ = {}.__delitem__
  else:
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

  def __repr__(self):
    if self.keys():
      L = list(self.keys())
      L = [str(L[i]) for i in range(len(L))]
      m = max(map(len, L)) + 1
      f = ''
      for k, v in self.items():
        if isinstance(v,dict):
          f = f + '\n'+str(k).rjust(m) + ': ' + repr(k) + ' class'
        else:
          unitStr=''
          if k in unit:
            unitStr = ' ['+unit[k]+']'
          f = f + '\n'+str(k).rjust(m) + ': ' + repr(v) + unitStr
      return f
    else:
      return self.__class__.__name__ + "()"

  def find_key(self,val):
    if val==None :
      return
    for k in self.keys():
      if self[k]==val:
        return k

unit = dictClass({'length'       :'m',
                  'mass'         :'GeV',
                  'charge'       :'e',
                  'current'      :'A',
                  'energy'       :'GeV',
                  'ekinetic'     :'GeV',
                  'fint'         :'1'
                 })

    
      
#%%============================================================================
#                                global settings                               
#============================================================================== 
class elemClass(dictClass):
  """
  command container
  """
  def update(self):
    """
    update member variables according to each command data structure
    """
  def str(self):
    f = self.name+': '+self.elem+', '
    for k, v in self.items():
      if k not in ['name','elem']:
        f = f + k +'='+str(v).lower()+', '
    return f[:-2]
  def _input_multipleOptionals(self,keys,vals):
    nvals=0
    for val in vals:
      if val!=None:
        nvals=nvals+1
    if nvals>1:
      raise ValueError('data conflict: only one of '+str(keys)[1:-1]+' must be present.')
    else:
      for i,val in enumerate(vals):
        if val!=None:
          self[keys[i]]=val
  def _check_multipleOptionals(self,keys):
    keyIn= []
    for key in keys:
      if key in self.keys():
        keyIn.append(key)
    if len(keyIn)>1:
      raise ValueError('data conflict: only one of '+str(keys)[1:-1]+' must be present. But '+str(keyIn)[2:-2]+' are present.')
    if len(keyIn)==0:
      UserWarning('data conflict: at least one of '+str(keys)[1:-1]+' must be present.')
          

class beam(elemClass):
  """
  beam parameter setting
  """
  def __init__(self,name='beam',particle=None,mass=None,charge=None,energy=None,ekinetic=None):
    self.elem  = 'beam'
    self.name = name
    self._input_multipleOptionals(['energy','ekinetic'],[energy,ekinetic])
    self._input_multipleOptionals(['particle','mass'],[particle,mass])
    self._input_multipleOptionals(['particle','charge'],[particle,charge])
  def update(self):
    self._check_multipleOptionals(['particle','mass'])
    self._check_multipleOptionals(['particle','charge'])
    self._check_multipleOptionals(['energy','ekinetic'])
  
  
class units(elemClass):
  """
  beam parameter setting
  """
  def __init__(self,name='units',type='static'):
    self.elem  = 'units'
    self.name = name
    self.type = type


class globaldefaults(elemClass):
  """
  beam parameter setting
  """
  def __init__(self,name='setdefaults',lfrngsbend=1,tfrngsbend=1,lfrngquad=0,tfrngquad=0,driftexact=0):
    self.elem = 'globaldefaults'
    self.name = name
    self.lfrngsbend = lfrngsbend
    self.tfrngsbend = tfrngsbend
    self.lfrngquad  = lfrngquad
    self.tfrngquad  = tfrngquad
    self.driftexact = driftexact

def elemList(beam):
  return [beam,units(),globaldefaults()]
    
#%%============================================================================
#                                  commands                               
#==============================================================================
class autotrack(elemClass):
  """
  particle tracking
  """
  def __init__(self,name='autotrack',type='symplectic',order='5'):
    self.elem = 'autotrack'
    self.name = name
    if type not in ['symplectic','taylor']:
      raise ValueError(type)
    self.type = type
    self.order = order
  def update(self):
    if self.type[:-1] in ['symplectic','taylor']:
      try:
        self.order=int(self.type[-1])
      except: 
        ValueError(self.type)
      self.type=self.type[:-1]
  def str(self):
    return self.name+':autotrack, type='+self.type+str(self.order)
    

class particledump(elemClass):
  """
  write particle data on file
  """
  def __init__(self,name='particledump',file='rays.out',precision=9,close=False,flush=True,nunits=0):
    self.elem = 'particledump'
    self.elem = 'particledump'
    self.name = name
    self.file = file
    self.precision = precision
    self.close = close
    self.flush = flush
    self.nunits = nunits 

    
class iden(elemClass):
  """
  identity map
  """
  def __init__(self,name='clear'):
    self.elem = 'iden'
    self.name = name
clear = iden()


class end(elemClass):
  """
  finish MLI
  """
  def __init__(self,name='fin'):
    self.elem = 'end'
    self.name = name
fin = end()


class raytrace(elemClass):
  """
  read particle data
  """
  def __init__(self,name='raytrace',file1='rays.in',type='readonly'):
    self.elem = 'raytrace'
    self.name = name
    self.file1 = file1
    self.type = type
  
class ptm(elemClass):
  """
  print map
  inputs:
    matrix = ?
    poly = polynomial order
  """
  def __init__(self,name='ptm',matrix=3,poly=3):
    self.elem = 'ptm'
    self.name = name
    self.matrix = matrix
    self.poly = poly

    
class aim(elemClass):
  """
  optimization target
  inputs:
    infile = target info file number
    job = 1:specify aim only, 2:specify aim and targets, 3:specify aim, target and weights
  """
  def __init__(self,name='aim',infile=10,job=2):
    self.elem = 'aim'
    self.name = name
    self.infile = infile
    self.job = job
    
    
class vary(elemClass):
  """
  optimization knob for adjust
  inputs:
    infile = knob info file number
    job = 1:NV = NA, 2:NV <= NA, 3 others  (where NV = number of knobs, NA = number of aims)
  """
  def __init__(self,name='vary',infile=10,job=2):
    self.elem = 'vary'
    self.name = name
    self.infile = infile
    self.job = job
    
  
class tasm(elemClass):
  """
  calculate (chromatic, anharmonic) tune 
  """
  def __init__(self,name='tasm',iopt=1,delta=0.0,idata=1,ipmaps=0,isend=3,iwmaps=0):
    self.elem = 'tasm'
    self.name = name
    self.iopt = iopt
    self.delta = delta
    self.idata = idata
    self.ipmaps = ipmaps
    self.isend = isend
    self.iwmaps = iwmaps
  
    
class stm(elemClass):
  """
  store map in MLI memory in location nmap
  """
  def __init__(self,name='stm',nmap=1):
    self.elem = 'stm'
    self.name = name
    self.nmap = nmap
  
  
class gtm(elemClass):
  """
  get map from MLI memory in location nmap
    inputs: 
    nmap = MLI memory location (int)
    iopt = 1:concatenate, 2=replace map
  """
  def __init__(self,name='gtm',nmap=1,iopt=2):
    self.elem = 'gtm'
    self.name = name
    self.nmap = nmap
    self.iopt = iopt
    
  
class tmo(elemClass):
  """
  output transfer map
    inputs: 
    ifile = file number s.t. fort.ifile
  """
  def __init__(self,name='tmo',ifile=17):
    self.elem = 'tmo'
    self.name = name
    self.ifile = ifile
#%%============================================================================
#                                   lattice                               
#==============================================================================
class line(elemClass):
  """
  build line from elements
  inputs: 
    elemList = list of elements
  """
  def __init__(self,name='lattice',elemList=[]):
    self.elem  = 'line'
    self.name = name
    #self.list = elemList.copy()
    self.list=[]
    for item in elemList:
      if isinstance(item, str):
        self.list.append(item)
      else:
        try:
          self.list.append(item.name)
        except TypeError:
          print('following is not MLI element')
          print(item)
          
  def str(self):
    f = self.name+', line=( '
    k=0
    for item in self.list:
      k=k+1
      f = f+ item + ' '
#       if isinstance(item, str):
#         f = f+ item + ' '
#       else:
#         try:
#           f = f+item.name + ' '
#         except TypeError:
#           print('following is not MLI element')
#           print(item)
      if k>=5:
        k=0
        f = f+ '& \n'
    f = f+')'
    return f
  
  
class tmi(elemClass):
  """
  get map from file
  inputs: 
    iopt = 1:concatenate, 2=replace map
    nopt = 1:rewind file
    iflie = read from fort.ifile 
    nskip = 
  """
  def __init__(self,name='tmi',ifile=16,iopt=2,nopt=1,nskip=0):
    self.elem = 'tmi'
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
  
  
class nlinsert(elemClass):
  """
  nonlinear insert element
  inputs: 
    iopt = 1:concatenate, 2=replace map
    nopt = 1:rewind file
    iflie = read from fort.ifile 
    nskip = 
  """
  def __init__(self,name='nlinsert',zstart=0.0,zend=1.8,steps=1000,zlen=1.8,k=1.45446332708327,tau=-0.4,c=0.01):
    self.elem  = 'nlinsert'
    self.name = name
    self.zstart = zstart
    self.zend = zend
    self.steps = steps
    self.zlen = zlen
    self.k = k
    self.c = c
    self.tau = tau
  
  
class monitor(elemClass):
  """
  monitor
  """
  def __init__(self,name='monitor'):
    self.elem  = 'monitor'
    self.name = name  
    

class marker(elemClass):
  """
  marker
  """
  def __init__(self,name='marker'):
    self.elem  = 'marker'
    self.name = name  
    

class drift(elemClass):
  """
  drift
  inputs:
    l  = length
  """
  def __init__(self,name='drift',l=0.1):
    self.elem  = 'drift'
    self.name = name
    self.l = l
    
    
class quadrupole(elemClass):
  """
  quadrupole
  inputs:
    l  = length
    k1 = gradient
  """
  def __init__(self,name='quadrupole',l=0.1,g1=None,k1=None):
    self.elem  = 'quadrupole'
    self.name = name
    self.l = l
    self._input_multipleOptionals(['g1','k1'],[g1,k1])
  def update(self):
    self._check_multipleOptionals(['g1','k1'])
      
      
class vkicker(elemClass):
  """
  vkicker
  inputs:
    l  = length
    kick = kick strengh
  """
  def __init__(self,name='vkicker',l=0.1,kick=0):
    self.elem  = 'vkicker'
    self.name = name
    self.l = l
    self.kick = kick
  
  
class sextupole(elemClass):
  """
  sextupole
  inputs:
    l  = length
    g2 = sextupole strengh
  """
  def __init__(self,name='sextupole',l=0.1,g2=None):
    self.elem  = 'sextupole'
    self.name = name
    self.l = l
    if g2!=None:
      self.g2 = g2
      
      
class thlm(elemClass):
  """
  thin low-order multiple
  inputs:
    k1l = integrated quadrupole strengh
    k2l = integrated sextupole strengh
    k3l = integrated octupole strengh
  """
  def __init__(self,name='thlm',k1l=None,k2l=None,k3l=None):
    self.elem  = 'thlm'
    self.name = name
    if k1l!=None:
      self.k1l = k1l
    if k2l!=None:
      self.k2l = k2l
    if k3l!=None:
      self.k3l = k3l
    
      
class sbend(elemClass):
  """
  sector dipole
  inputs:
    l = length
    angle 
    hgap
    fint = fringe field integration parameter
  """
  def __init__(self,name='sbend',l=0.1,angle=0,hgap=0,fint=0):
    self.elem  = 'sbend'
    self.name = name
    self.l = l
    self.angle = angle
    self.hgap = hgap
    self.fint = fint
      
      
class dipedge(elemClass):
  """
  dipole edge
  inputs:
    e1 = edge angle
    h 
    hgap
    fint = fringe field integration parameter
  """
  def __init__(self,name='dipedge',e1=0,h=0,hgap=0,fint=0):
    self.elem  = 'dipedge'
    self.name = name
    self.e1 = e1
    self.h = h
    self.hgap = hgap
    self.fint = fint
  
  
class srfc(elemClass):
  """
  thin rf cavity
  """
  def __init__(self,name='srfc',volts=0,freq=0):
    self.elem  = 'srfc'
    self.name = name
    self.volts = volts
    self.freq = freq