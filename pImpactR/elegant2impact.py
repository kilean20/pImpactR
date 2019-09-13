import os
import numpy as np
from copy import deepcopy as copy
import re

__all__ = ['elegant2impact']


def _defaultKeyVal(d,k,v):
  if k in d.keys():
    return d[k]
  else:
    return v
    

class _dictClass(dict):
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

      
unit = _dictClass({'length'       :'m',
                  'mass'         :'eV',
                  'charge'       :'e',
                  'current'      :'A'
                 })   
      
####################################
#   read elegant lattice file
####################################
def _isint(s):
  try:
    int(s)
    return True
  except ValueError:
    return False

  
def _isfloat(s):
  try:
    float(s.replace('D','E',1).replace('d','e',1))
    return True
  except ValueError:
    return False
  
def _isvariable(s,variables):
  if s in variables.keys():
    return True
  else:
    return False
  
def _str2val(s,variables):
  if _isint(s):
    return int(s)
  elif _isfloat(s):
    return float(s.replace('D','E',1).replace('d','e',1))
  elif _isvariable(s,variables):
    return variables[s]
  else:
    return s

def _file2rawlines(fname):
  lines = []
  tmp = ''
  variables = {}
  with open(fname) as f:
    flagLine = False
    for line in f:
      line=line.replace('\n','')
      if line.lstrip().find('!')==0 or line.lstrip()=='' or line.lstrip()[:3]=='>>>':
        continue
      i = line.find('!')
      if i!=-1:
        line = line[:i]  
      if ('line' in line or 'LINE' in line)  and '(' in line:
        flagLine = True
      if flagLine:
        line=line.replace('&','')
        i=line.find(')')
        if i!=-1:
          flagLine = False
          tmp = tmp + line
          lines.append(copy(tmp))
          tmp = ''
        else:
          tmp = tmp + line
      else:
        line=line.replace(' ','')
        lines.append(line)
  return lines

def _rawlines2var(raw):
  var = {}
  raw2 = raw.copy()
  k=0
  for i,line in enumerate(raw):
    icomma = line.find(',')
    icolon = line.find(':')
    iequal = line.find('=')
    if icomma==-1 and icolon==-1 and iequal!=-1:
      var[line[:iequal]]=_str2val(line[iequal+1:],var)
      del raw2[i-k]
      k=k+1
  return var,raw2
  
  
def _rawlines2elem(raw,variables):
  elems=[]
  raw2 = raw.copy()
  k=0
  for i,line in enumerate(raw):
    icolon = line.find(':')
    if icolon!=-1:
      del raw2[i-k]
      k=k+1
      istart= line.find('{')
      if istart==-1:
        line = re.split(':|,|=',line)
      else:
        iend  = line.find('}')
        line = re.split(':|,|=',line[:istart]) + [line[istart:iend+1]] + re.split(':|,|=',line[iend+1:])
        line = list(filter(('').__ne__, line))
      name = line[0]
      elem  = line[1]
      f = _dictClass()
      f.elem = elem
      f.name = name
      if not elem in ['drif','kquad','ksext','csbend']:
        print(elem + ' is not recognized. skipping...')
        continue
      if len(line)>2:
        for i in range(2,len(line),2):
          f[line[i]]=_str2val(line[i+1],variables)
      f.update()
      elems.append(f)
  return elems,raw2


class _line(_dictClass):
  """
  build line from elements
  inputs: 
    elemList = list of elements
  """
  def __init__(self,name='lattice',elemList=[]):
    self.elem  = 'LINE'
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
          print('following is not elegant element')
          
  def str(self):
    f = self.name+', line=( '
    k=0
    for item in self.list:
      k=k+1
      f = f+ item + ' '
      if k>=5:
        k=0
        f = f+ '& \n'
    f = f+')'
    return f

def _rawlines2lattice(raw,elems):
  raw2=[]
  for line in raw:
    line=line.replace('&',' ')
    line=line.replace('(',' ')
    line=line.replace(')',' ')
    line=line.replace('=',' ')
    line=line.replace(',',' ')
    raw2 = raw2 + line.split()
  iElems = [i for i in range(len(raw2))]
  ilattices = []
  elemList = raw2.copy()
  k=0
  for i,item in enumerate(raw2):
    if item in ['line','LINE']:
      iElems.remove(i)
      iElems.remove(i-1)
      ilattices.append(i)

  for i in iElems:
    for elem in elems:
      if elem.name == raw2[i]:
        raw2[i]=elem

  nlattices = len(ilattices)
  lattices = [0]*nlattices
  i=-1
  for i in range(nlattices-1):
    lattices[i]=_line(name=raw2[ilattices[i]-1],elemList=raw2[ilattices[i]+1:ilattices[i+1]-1])
  lattices[i+1]=_line(name=raw2[ilattices[i+1]-1],elemList=raw2[ilattices[i+1]+1:])
  return lattices

def _rawlines_to_elemList_n_latticeList(rawlines):
  var,raw2=_rawlines2var(rawlines)
  elemList,raw2= _rawlines2elem(raw2,var)
  latticeList  = _rawlines2lattice(raw2,elemList)
  return elemList,latticeList

def _readElegant(fname='ele.lat'):
  rawlines = _file2rawlines(fname)
  elemList,latticeList = _rawlines_to_elemList_n_latticeList(rawlines)
  return elemList,latticeList


###############################################
#   convert elegant lattice dict to impact
###############################################


def elegant2impact(fname,mass,kinetic_energy,charge,min_sext_L=0.02):
  from data import beam as _beam
  from impactIO import getElem
  beam=_beam()
  elemList,elegLine = _readElegant(fname=fname)
  beam.distribution.distribution_type = 'ReadFile'
  beam.mass = mass
  beam.kinetic_energy = kinetic_energy
  beam.charge = charge
  beam.multi_charge.q_m[0] = beam.charge/beam.mass
  
  cLight = 299792458.0
  rel_gam = kinetic_energy/mass+1.0
  rel_bet = np.sqrt((rel_gam+1.0)*(rel_gam-1.0))/rel_gam
  rigidity = cLight/(mass*rel_bet*rel_gam)
  
  
  elegLine = elegLine[0]
  fList = []
  for name in elegLine.list:
    f=None
    for item in elemList:
      if item.name==name:
        elem=item.elem
        if   elem == 'csbend':
          f = getElem('dipole')
          f.file_id = 350
          f.length = item.L
          if 'HGAP' not in item.keys():
            f.pipe_radius = 1.0
          else:
            f.pipe_radius = 2.0*item.HGAP
          if 'FINT' not in item.keys():
            f.fringe_field_integration = 0.0
          else:
            f.fringe_field_integration = item.FINT
          f.bending_angle = item.angle
          if 'e1' not in item.keys():
            f.entrance_angle = 0.0
          else:
            f.entrance_angle = item.e1
          if 'e2' not in item.keys():
            f.exit_angle = 0.0
          else:
            f.exit_angle = item.e2
          f=[f]
        elif elem == 'drif':
          f = getElem('drift')
          f.length = item.L
          f=[f]
        elif elem == 'kquad':
          f = getElem('quad')
          f.length = item.L
          f.B1 = item.K1#/rigidity
          f=[f]
        elif elem == 'ksext':
          L = item.L
          N = int(np.ceil(L/min_sext_L))
          print('elegant ksext of length '+str(L)+' found. Converting to '+str(N)+' drift-kick-drift thin element')
          kick = getElem('multipole_thin')
          d0 = getElem('drift')
          d1 = getElem('drift')
          kick.KL_sext = item.K2*L/N
          d0.length = 0.5*L/N
          d1.length = L/N
          f=[d0,kick]
          for i in range(N-1):
            f = f + [copy(d1),copy(kick)]
          f=f+[copy(d0)]              
          
        else:
          print(elem + ' is not recognized. skipping...')
        break
    if f!= None:
      fList=fList+f
    for elem in fList:
      if 'length' in elem:
        elem.n_sckick = int(np.ceil(elem.length*20))
  return beam, fList
  
  