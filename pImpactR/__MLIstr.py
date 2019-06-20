import numpy as np
import os


# def _pad_equal(s):
#   flagRedo = True
#   while flagRedo:
#     flagRedo = False
#     for i,item in enumerate(s):
#       if item=='=':
#         if s[i-1]==' ' and s[i+1]!=' ':
#           s = s[:i-1]+'= '+s[i+1:]
#           flagRedo = True
#           break
#         elif s[i-1]==' ' and s[i+1]==' ':
#           s = s[:i-1]+'='+s[i+1:]
#           flagRedo = True
#           break
#         elif s[i+1]!=' ':
#           s = s[:i]+'= '+s[i+1:]
#           flagRedo = True
#           break
#   return s


# def _DtoE(word):
#   if 'D' in word or 'd' in word: 
#     try:
#       temp = float(word.replace('D','E',1).replace('d','e',1))
#       return str(temp)
#     except:
#       return word
#   else:
#     return word

  
# def _readraw(fname):
#   try:
#     fin = open(fname,'r')
#     dataList  = fin.readlines()
#     fin.close()
#   except:
#     print(( "  ERROR! Can't open file '" + fname + "'"))
#     return False

#   i=0
#   while i < len(dataList):
#     if dataList[i].lstrip()=='' or dataList[i].lstrip().startswith('!'):
#       del dataList[i]
#       i=i-1
#     else:
#       index = dataList[i].lstrip().find('!')
#       if index==-1:
#           dataList[i]=dataList[i].strip()+'\n'
#       else:
#           dataList[i]=dataList[i].lstrip()[:index].rstrip()+'\n'
#     i=i+1
#   for line in dataList:
#   dataList  = [_pad_equal(line).split() for line in dataList ]
#   for i in range(0,len(dataList)):
#     for j in range(0,len(dataList[i])):
#       if dataList[i][j] == '/':
#         del dataList[i][j:]
#     for j in range(0,len(dataList[i])):
#       dataList[i][j] = _DtoE(dataList[i][j])
#   return dataList


def getMenu(fname=None,particle='electron',energy=0.150511006,fringfield=False,driftexact=False):
  if fname==None:
    f = 'beam: beam, particle='+particle+',energy='+str(energy)+' \n'
    f = f+'units, type=static \n'
    if fringfield:
      f = f+'setdefaults: globaldefaults, lfrngsbend=1,tfrngsbend=1,lfrngquad=1,tfrngquad=1,driftexact='+str(int(driftexact))
    else:
      f = f+'setdefaults: globaldefaults, lfrngsbend=0,tfrngsbend=0,lfrngquad=0,tfrngquad=0,driftexact='+str(int(driftexact))
    return f+' \n'
  else:
    with open(fname) as f:
      return f.read() + '\n\n'
          
    
def getLabor(works):
  f = '#labor \n'
  for item in works:
    f = f + item + '\n'
  return f


def run(nCore=None):
  if nCore==None:
    os.system('mli.x')
  else:
    os.system('mpirun -n '+str(nCore)+' mli.x')
  