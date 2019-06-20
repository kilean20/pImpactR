import os
import numpy as np
from data import dictClass
import pandas as pd
import MLI_getElem as getElem
from copy import deepcopy as copy

def run(nCore=None):
  if nCore==None:
    os.system('x')
  else:
    os.system('mpirun -n '+str(nCore)+' x')

def cmd2str(commandList):
  f = '\n'
  for item in commandList:
    if isinstance(item, str):
      f = f+ item + ' \n'
    else:
      try:
        f = f+item.str() + ' \n'
      except TypeError:
        print('following is not MLI command')
        print(item)
  return f

  
class buildMenu(dictClass):
  def __init__(self,commandList=[]):
    self.list = commandList
  def str(self):
    return '#menu \n' + cmd2str(self.list)

  
class buildLabor(dictClass):
  def __init__(self,commandList=[]):
    self.list = commandList
  def str(self):
    f = '\n#labor \n'
    for item in self.list:
      if isinstance(item, str):
        f = f+ item + ' \n'
      elif isinstance(item, int):
        f = f+ str(item) + '*'
      else:
        try:
          f = f+item.name + ' \n'
        except TypeError:
          print('following is not MLI command')
          print(item)
    return f

def writeInputfile(menu,labor,fname='in'):
  with open(fname,'w') as f:
    f.write(menu.str())
    f.write(labor.str())
  
#%%============================================================================
#                                   read outputs                               
#==============================================================================
def readMatrix(fname='out'):
  M = np.zeros([6,6])
  flagIn = False
  with open(fname) as f:
    for i,line in enumerate(f):
      if flagIn:
        if line=='\n':
          break
        items = line.strip().split()
        M[int(items[0])-1,int(items[1])-1]=float(items[2])
      if 'nonzero matrix elements in full precision:' in line:
        flagIn = True      
  return pd.DataFrame(M,index=range(1,7),columns=range(1,7))

def readGeneratingPolynomial(fname='out'):
  flagIn = False
  flagIn2 = False
  Ind1 = []
  Ind2 = []
  Val  = []
  with open(fname) as f:
    for i,line in enumerate(f):
      if flagIn2:
        if line=='\n':
          break
        #Ind1.append(line[2:22])
        Ind1.append(int(line[4:7]))
        Ind2.append(line[9:22])
        Val.append(float(line[23:-1]))
        
      if flagIn:
        if line=='\n':
          flagIn2 = True
      if 'nonzero elements in generating polynomial are :' in line:
        flagIn = True
  return pd.DataFrame({'exponents':Ind2,'value':Val},index=Ind1)

def getTBT(npt,nturn,fname='rays.out'):
    TBT = np.loadtxt(fname)
    dim = TBT.ndim
    if dim==1:
      return TBT[1:-1]
    else:
      dummy,nCol = TBT.shape
      if nCol>7:
        TBT = TBT[:npt*nturn,1:-1]
      else:
        TBT = TBT[:npt*nturn,:-1]
      out = np.zeros([npt,nturn,6])
      for i in range(nturn):
        out[:,i,:] = TBT[i*npt:(i+1)*npt,:].reshape([npt,6])
      return out

    
#%%============================================================================
#                                   read inputs                               
#==============================================================================
def _isint(value):
  try:
    int(value)
    return True
  except ValueError:
    return False

  
def _isfloat(value):
  try:
    float(value.replace('D','E',1).replace('d','e',1))
    return True
  except ValueError:
    return False

  
def _isvariable(value,variables):
  if value in variables.keys():
    return True
  else:
    return False

  
def _str2val(s,variables):
  if _isint(s):
    return int(s)
  elif _isfloat(s):
    return float(s)
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
      if 'line' in line:
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
        try:
          i=line.index('!')
        except:
          i=len(line)
        if i>0:
          lines.append(line[:i])
  return lines


def _rawlines2menu(lines):
  variables = {}
  menu = []
  for line in lines:
    try:
      i = line.index(':')
      name = line[:i]
      line=line[i+1:]
      i = line.find(',')
      if i>0:
        Elemtype = line[:i]
        line=line[i+1:]
      else:
        Elemtype = line[:]
      tmp = {}
      nkey = 0
      for i in range(len(line)):
        if line[i]=='=':
          nkey=nkey+1
      for ikey in range(nkey):
        i = line.index('=')
        if ikey<=nkey-1:
          j = line.index(',')
          tmp[line[:i]]=line[i+1:j]
          line = line[j+1:]
        else:
          tmp[line[:i]]=line[i+1:]
      if Elemtype == 'sbend':
        f = getElem.sbend(name=name)
      elif Elemtype == 'beam':
        f = getElem.beam(name=name)
      elif Elemtype == 'units':
        f = getElem.units(name=name)
      elif Elemtype == 'globaldefaults':
        f = getElem.globaldefaults(name=name)
      elif Elemtype == 'autotrack':
        f = getElem.autotrack(name=name)
      elif Elemtype == 'particledump':
        f = getElem.particledump(name=name)
      elif Elemtype == 'iden':
        f = getElem.iden(name=name)
      elif Elemtype == 'end':
        f = getElem.end(name=name)
      elif Elemtype == 'raytrace':
        f = getElem.raytrace(name=name)
      elif Elemtype == 'ptm':
        f = getElem.ptm(name=name)
      elif Elemtype == 'tasm':
        f = getElem.tasm(name=name)
      elif Elemtype == 'stm':
        f = getElem.stm(name=name)
      elif Elemtype == 'gtm':
        f = getElem.gtm(name=name)
      elif Elemtype == 'tmo':
        f = getElem.tmo(name=name)
      elif Elemtype == 'tmi':
        f = getElem.tmi(name=name)
      elif Elemtype == 'nlinsert':
        f = getElem.nlinsert(name=name)
      elif Elemtype == 'monitor':
        f = getElem.monitor(name=name)
      elif Elemtype == 'marker':
        f = getElem.marker(name=name)
      elif Elemtype == 'drift':
        f = getElem.drift(name=name)
      elif Elemtype == 'quadrupole':
        f = getElem.quadrupole(name=name)
      elif Elemtype == 'vkicker':
        f = getElem.vkicker(name=name)
      elif Elemtype == 'sextupole':
        f = getElem.sextupole(name=name)
      elif Elemtype == 'sbend':
        f = getElem.sbend(name=name)
      elif Elemtype == 'dipedge':
        f = getElem.dipedge(name=name)
      elif Elemtype == 'srfc':
        f = getElem.srfc(name=name)
      for k,v in tmp.items():
        f[k] = _str2val(v,variables)
      print(f.items())
    except:
      i = line.index('=')
      v = line[i+1:]
      variables[line[:i]] = _str2val(v,variables)
      f=False
    if f:
      menu.append(f)
  return menu


def readInputfile(fname='in'):
  labor = []
  lines = _file2rawlines(fname)
  i_menu = -1
  for i_labor,line in enumerate(lines):
    if line[:6]=='#menu':
      i_menu = copy(i_labor)
      break
    if line[:6]=='#labor':
      break
  menu = _rawlines2menu(lines[i_menu+1:i_labor])
  return menu


  