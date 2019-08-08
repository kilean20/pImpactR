import os
import numpy as np
from data import dictClass
import pandas as pd
import MLI_getElem as getElem
from copy import deepcopy as copy
import re

def run(nCore=None):
  if nCore==None:
    os.system('time mli.x > mli.log')
  else:
    os.system('time mpirun -n '+str(nCore)+' mli.x > mli.log')

def elem2str(elemList):
  f = '\n'
  for item in elemList:
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
  def __init__(self,elemList=[],latticeList=[]):
    self.elemList = elemList
    self.latticeList = latticeList
  def str(self):
    return '#menu \n' + elem2str(self.elemList) + elem2str(self.latticeList)


class buildLabor(dictClass):
  def __init__(self,elemList=[]):
    self.list = elemList
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
  
#%%============================================================================
#                                   read outputs                               
#==============================================================================
def readTransferMap(fname='mli.out'):
  M = np.zeros([6,6])
  Ind1 = []
  Ind2 = []
  Val  = []
  iM = -1
  iG = -1
  with open(fname) as f:
    lines = f.readlines()
  if fname[:4]=='fort':
    flagM = False
    flagG = False
    for i,line in enumerate(lines):
      line = line.strip().split()
      lines[i] = line
      if len(line) == 3 and not flagM:
        iM = i
        flagM = True
        flagG = False
      if len(line) == 2 and not flagG:
        iG = i
        flagM = False
        flagG = True
    # --- get Matrix --- 
    if iM != -1:
      for line in lines[iM:iG]:
        M[int(line[0])-1,int(line[1])-1]=float(line[2])
    M = pd.DataFrame(M,index=range(1,7),columns=range(1,7))
    # --- get generating polynomial ---
    if iG != -1:
      for line in lines[iG:]:
        Ind1.append(int(line[0]))
        Val.append(float(line[1]))
    G = pd.DataFrame({'GP':Val},index=Ind1)
  elif fname=='mli.out':
    for i,line in enumerate(lines):
      if 'nonzero matrix elements in full precision:' in line:
        iM = i
      if 'nonzero elements in generating polynomial are :' in line:
        iG = i
    # --- get Matrix --- 
    if iM != -1:
      for line in lines[iM+1:]:
        if line=='\n':
            break
        items = line.strip().split()
        M[int(items[0])-1,int(items[1])-1]=float(items[2])
    M = pd.DataFrame(M,index=range(1,7),columns=range(1,7))
    # --- get generating polynomial ---
    if iG != -1:
      with open(fname) as f:
        for line in lines[iG+2:]:
          if line=='\n':
            break
          #Ind1.append(line[2:22])
          Ind1.append(int(line[4:7]))
          Ind2.append(line[9:22])
          Val.append(float(line[23:-1]))
    G = pd.DataFrame({'exponents':Ind2,'GP':Val},index=Ind1)
  else:
    raise ValueError('fname')
  return M,G  

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
#                                   write inputs                               
#==============================================================================
def writeInputfile(elemList,latticeList,labor,fname='mli.in'):
  menu = buildMenu(elemList,latticeList)
  with open(fname,'w') as f:
    f.write(menu.str())
    try:
      f.write(labor.str())
    except:
      laborStr='#labor\n'
      for item in labor:
        laborStr=laborStr+item+'\n'
      f.write(laborStr)
    
    
#%%============================================================================
#                                   read inputs                               
#==============================================================================
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
      if 'line' in line and '(' in line:
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
      if elem   == 'beam':
        f = getElem.beam(name=name)
      elif elem == 'units':
        f = getElem.units(name=name)
      elif elem == 'globaldefaults':
        f = getElem.globaldefaults(name=name)
      elif elem == 'autotrack':
        f = getElem.autotrack(name=name)
      elif elem == 'particledump':
        f = getElem.particledump(name=name)
      elif elem == 'iden':
        f = getElem.iden(name=name)
      elif elem == 'end':
        f = getElem.end(name=name)
      elif elem == 'raytrace':
        f = getElem.raytrace(name=name)
      elif elem == 'ptm':
        f = getElem.ptm(name=name)
      elif elem == 'tasm':
        f = getElem.tasm(name=name)
      elif elem == 'stm':
        f = getElem.stm(name=name)
      elif elem == 'gtm':
        f = getElem.gtm(name=name)
      elif elem == 'tmo':
        f = getElem.tmo(name=name)
      elif elem == 'tmi':
        f = getElem.tmi(name=name)
      elif elem == 'vary':
        f = getElem.vary(name=name)
      elif elem == 'aim':
        f = getElem.aim(name=name)
      elif elem == 'nlinsert':
        f = getElem.nlinsert(name=name)
      elif elem == 'monitor':
        f = getElem.monitor(name=name)
      elif elem == 'marker':
        f = getElem.marker(name=name)
      elif elem == 'drift':
        f = getElem.drift(name=name)
      elif elem == 'quadrupole':
        f = getElem.quadrupole(name=name)
      elif elem == 'vkicker':
        f = getElem.vkicker(name=name)
      elif elem == 'sextupole':
        f = getElem.sextupole(name=name)
      elif elem == 'thlm':
        f = getElem.thlm(name=name)
      elif elem == 'sbend':
        f = getElem.sbend(name=name)
      elif elem == 'dipedge':
        f = getElem.dipedge(name=name)
      elif elem == 'srfc':
        f = getElem.srfc(name=name)
      else:
        print(elem + ' is not recognized. skipping...')
        continue
      if len(line)>2:
        for i in range(2,len(line),2):
          f[line[i]]=_str2val(line[i+1],variables)
      f.update()
      elems.append(f)
  return elems,raw2


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
    if item == 'line':
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
    lattices[i]=getElem.line(name=raw2[ilattices[i]-1],elemList=raw2[ilattices[i]+1:ilattices[i+1]-1])
  lattices[i+1]=getElem.line(name=raw2[ilattices[i+1]-1],elemList=raw2[ilattices[i+1]+1:])
  return lattices
    

def _rawlines_to_elemList_n_latticeList(rawlines):
  var,raw2=_rawlines2var(rawlines)
  elemList,raw2=_rawlines2elem(raw2,var)
  latticeList  =_rawlines2lattice(raw2,elemList)
  return elemList,latticeList


def readInputfile(fname='mli.in'):
  labor = []
  rawlines = _file2rawlines(fname)
  i_menu = -1
  for i_labor,line in enumerate(rawlines):
    if line=='#menu':
      i_menu = copy(i_labor)
    if line=='#labor':
      break
  elemList,latticeList = _rawlines_to_elemList_n_latticeList(rawlines[i_menu+1:i_labor])
  labor= rawlines[i_labor+1:]
  return elemList,latticeList,labor
    
    
#%%============================================================================
#                             Lattice Menipulate                    
#==============================================================================
def sext2thin(elemList,latticeList,brho=None):
  newList=[]
  sextupoles={}
  for item in elemList:
    if item.elem == 'sextupole':
      thin_drift = getElem.drift(name='drift_thin_'+item.name,l=0.5*item.l)
      if 'k2' in item:
        thin_multi = getElem.thlm(name=item.name,k2l=item.l*item.k2)
      elif 'g2' in item:
        if brho!=None:
          thin_multi = getElem.thlm(name=item.name,k2l=item.l*item.g2*2.0/brho)
        else:
          raise ValueError('argument brho is needed to convert g2 to k2')
      else:
        UserWarning('the sextupole '+item.name+' does not have g2 or k2 defined.')
        thin_multi = getElem.thlm(name=item.name,k2l=0.0)
      tmp = [thin_drift,thin_multi].copy()
      sextupoles[item.name]=[thin_drift.name,item.name,thin_drift.name]
    else:
      tmp = [item].copy()
    newList = newList+tmp
  newLatticeList = latticeList.copy()
  for lattice in newLatticeList:
    tmpLattice = lattice.list.copy()
    iSext = 0
    for i,item in enumerate(tmpLattice):
      if item in sextupoles.keys():
        lattice.list[i+iSext]=sextupoles[item][0]
        lattice.list.insert(i+iSext+1,sextupoles[item][1])
        lattice.list.insert(i+iSext+2,sextupoles[item][2])
        iSext=iSext+2
  return newList,newLatticeList
        
  
def removeElems(elem,elemList,latticeList):
  newList=[]
  names=[]
  for item in elemList:
    if item.elem == elem:
      names.append(item.name)
    else:
      newList.append(item)
  newLatticeList = latticeList.copy()
  for lattice in newLatticeList:
    tmpLattice = lattice.list.copy()
    iElem = 0
    for i,item in enumerate(tmpLattice):
      if item in names:
        lattice.list.pop(i-iElem)
        iElem=iElem+1
  return newList,newLatticeList
      
     