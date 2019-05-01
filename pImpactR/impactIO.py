from __future__ import print_function
from copy import deepcopy as copy
import numpy as np
import os
import data
import numbers
#====f2py routines=====
import read_phasespace as _read_pdata
import readTBT as _readTBT
#======================

#=======================================================================
#========= read turn by turn data ======================================
#=======================================================================
def intStr(eStr):
  return int(float(eStr))

def readTBT(fID, ke, mass, freq, nturn=None, path='.'):
  cwd=os.getcwd()
  os.chdir(path)
  if not os.path.isfile('TBT.'+str(fID)):
    print('can not find <TBT.'+str(fID)+'> file')
    return 'file error'
  
  if nturn==None:
    nturn,npt = _readTBT.get_tbtsize(fID)
  else:
    npt = _readTBT.get_tbtsize_npt(fID,nturn)
  tbt = _readTBT.get_tbtdata(fID,nturn,npt,ke,mass,freq)
  os.chdir(cwd)
  return tbt

def readTBTraw(fID, ke, mass, freq, nturn=None, path='.'):
  cwd=os.getcwd()
  os.chdir(path)
  
  if nturn==None:
    nturn,npt = _readTBT.get_rawtbtsize(fID)
  else:
    npt = _readTBT.get_rawtbtsize_npt(fID,nturn)
    
  tbt = _readTBT.get_rawtbtdata(fID,nturn,npt,ke,mass,freq)
  os.chdir(cwd)
  return tbt
  
def readTBT_integral(fID, nturn=None, path='.'):
  cwd=os.getcwd()
  os.chdir(path)
  if not os.path.isfile('TBT.integral.'+str(fID)):
    print('can not find <TBT.integral.'+str(fID)+'> file')
    return 'file error'
  
  if nturn==None:
    nturn,npt = _readTBT.get_tbtsize_integral(fID)
  else:
    npt = _readTBT.get_tbtsize_npt_integral(fID,nturn)

  tbt = _readTBT.get_tbtdata_integral(fID,nturn,npt)
  os.chdir(cwd)
  return tbt

def readTBT_integral_onMomentum(fID, nturn=None, path='.'):
  cwd=os.getcwd()
  os.chdir(path)
  
  if nturn==None:
    nturn,npt = _readTBT.get_tbtsize_integral_onmomentum(fID)
  else:
    npt = _readTBT.get_tbtsize_npt_integral_onmomentum(fID,nturn)

  tbt = _readTBT.get_tbtdata_integral_onmomentum(fID,nturn,npt)
  os.chdir(cwd)
  return tbt

#=======================================================================
#=======================================================================
#=======================================================================
def run(nCore=None,execfile='xmain'):
    """
    ierr = run(nCore=None,execfile='xmain')
    or
    ierr = run(beam,execfile='xmain')    
    run IMPACTz
    infer the source code directly and modify 
    to change executable path
    """  
    if nCore == None:
        return os.system(execfile+' > log.impact_std')
    elif isinstance(nCore,numbers.Integral) :
        return os.system('mpirun -n '+str(nCore) + ' ' + execfile +' > log.impact_std')
    else:
        return os.system('mpirun -n '+str(nCore.nCore_y*nCore.nCore_z) + ' ' + execfile +' > log.impact_std')
###############################################################################
###############################################################################
###                      IMPACT INPUT GENERATOR                             ###
###############################################################################
###############################################################################
#%%============================================================================
#                                    beam                               
#==============================================================================
def getBeam(dist_type='Waterbag') :
  beam = data.beam(dist_type)
  return beam

#%%============================================================================
#                                   lattice                                  
#==============================================================================
def getLattice() : 
  """
  f = getLattice()
  get a template (FODO lattce) of an lattice.  
  output 
      f = (list) list of elements
  """
  lattice = [getElem('loop'),getElem('quad'),getElem('drift'),getElem('quad'),getElem('drift')]
  lattice[3].B1 = -lattice[3].B1
  return lattice
#%%#================================element====================================
def getElem(type) : 
  """
  f = getElem(type)
  
  get a template of an element dictionary.  
  inquire data.elem_type for available types
  input 
      type = (str) element type. 
      
  output 
      f = (dict) element dictionary
  """
  if type not in list(data.elem_type.values()):
    raise AttributeError(type+'is not available type. See data.elem_type')
  elem = data.dictClass({'type':type})
  if type=='drift' :
    elem.length = 0.1
    elem.n_sckick = 1
    elem.n_map = 1
    elem.pipe_radius = 1.0
  elif type=='quad' :
    elem.length = 0.1
    elem.n_sckick = 1
    elem.n_map = 1
    elem.B1 = 10.0
    elem.file_id = 0
    elem.pipe_radius = 1.0
    elem.misalign_x = 0.0
    elem.misalign_y = 0.0
    elem.rotation_x = 0.0
    elem.rotation_y = 0.0
    elem.rotation_z = 0.0
  elif type=='const_focusing' :
    elem.length = 0.1
    elem.n_sckick = 1
    elem.n_map = 1
    elem.kx2 = 0.33333
    elem.ky2 = 0.33333
    elem.kz2 = 0.33333
    elem.pipe_radius = 1.0
  elif type=='solenoid' :
    elem.length = 0.1
    elem.n_sckick = 1
    elem.n_map = 1
    elem.Bz = 0.0
    elem.file_id = 0
    elem.pipe_radius = 1.0
    elem.misalign_x = 0.0
    elem.misalign_y = 0.0
    elem.rotation_x = 0.0
    elem.rotation_y = 0.0
    elem.rotation_z = 0.0
  elif type=='dipole' :
    elem.length = 0.1
    elem.n_sckick = 1
    elem.n_map = 1
    elem.bending_angle = 0.0
    elem.k1 = 0.0
    elem.file_id = 150
    elem.pipe_radius = 1.0
    elem.entrance_angle = 0.0
    elem.exit_angle = 0.0
    elem.entrance_curvature = 0.0
    elem.exit_curvature = 0.0
    elem.fringe_field_integration = 0.5
  elif type=='multipole_thin' :
    elem.n_sckick = 1
    elem.n_map = 1
    elem.KL_dipole = 0.0
    elem.KL_quad   = 0.0
    elem.KL_sext   = 0.0
    elem.KL_oct    = 0.0
    elem.KL_deca   = 0.0
    elem.KL_dodeca = 0.0
  elif type=='linear_matrix_map' :
    elem.nonlinear_insert_length = 1.0
    elem.nonlinear_insert_tuneAdvance = 0.3
    elem.tune_advance = 0.0
  elif type in ['nonlinear_insert','nonlinear_insert_smooth_focusing']  :
    elem.length = 1.0
    elem.n_sckick = 1
    elem.n_map = 1
    elem.strength_t = 0.0
    elem.transverse_scale_c = 0.003
    elem.pipe_radius = 1.0
    if type == 'nonlinear_insert':
      elem.tune_advance = 0.3
    else:
      elem.betx = 1.5
  elif type in ['TBT_integral','TBT_integral_onMomentum'] :
    elem.strength_t = 0.0
    elem.transverse_scale_c = 0.003
    elem.betx = 1.0
    elem.alfx = 0.0
    elem.file_id = 1000
    elem.pID_begin = 1
    elem.pID_end = 100
  elif type == 'TBT':
    elem.file_id = 1000
    elem.pID_begin = 1
    elem.pID_end = 100
  elif type == 'write_raw_ptcl':
    elem.file_id = 1000
    elem.format_id = 0
    elem.turn = 1
    elem.sample_period = 1
  elif type == 'pipe_override':
    elem.pipe_shape = 'rectangular'
    elem.xmax = 1.0
    elem.ymax = 1.0
  elif type=='loop':
    elem.turns = 1
  return elem

#%%============================================================================
#                                 input file
#==============================================================================
def readInputFile(fname='test.in'):
  """
  beam, lattice = readIMPACT(fname='test.in')
  read a IMPACT input file 
  output :
    beam : beam class containing beam and simulation control info
    lattice(list) element dictionaries
  """
  print('reading ImpactZ input file ('
        +data.bcolors.GREEN +fname 
        +data.bcolors.END+')' )
  raw = _readraw(fname)
  if raw == False:
      return
  beam = _str2beam(raw)
  """Lattice"""
  i=11
  print('  : lattice info ........................',end='') 
  lattice = []
  k=0
  for j in range(i,len(raw)):
    lattice.append(_str2elem(raw[j]))
    k+=1
  print('......done')
  return beam, lattice

def writeInputFile(beam,lattice,fname='test.in'):
  """
  write a IMPACT input file
  input 
      fname = (str) IMPACT input filename
      beam = (dict) beam dictionary
      lattce = (list) list of element dictionaries
  outfile
      fname
  """
  if sum(beam.multi_charge.n_particles) != beam.n_particles:
    #print('input error <- sum(beam.multi_charge.n_particles) not qual to beam.n_particles')
    if beam.multi_charge.n_states == 1:
      #print('  ... enforcing  beam.multi_charge.n_particles[0] to beam.n_particles')
      beam.multi_charge.n_particles[0]=beam.n_particles
    else:
      raise ValueError('program terminating...')
      
  if beam.multi_charge.n_states == 1 and beam.multi_charge.current[0] != beam.current :
    #print('input error <- beam.multi_charge.current[0] not qual to beam.current')
    #print('  ... enforcing  beam.multi_charge.current[0] to beam.current')
    beam.multi_charge.current[0] = beam.current
    
  beamStr = _beam2str(beam)
  for i in range(len(beamStr)):
    beamStr[i].append('\n')
    beamStr[i] = " ".join(beamStr[i])
    
  latticeStr = []
  for i in range(len(lattice)):
    latticeStr.append(_elem2str(lattice[i]))
    latticeStr[i].append('/')
    latticeStr[i].append('\n')
    latticeStr[i] = " ".join(latticeStr[i])
    
  f=open(fname,'w') 
  f.writelines(['!================= Beam & Control Parameters ================= \n'])
  f.writelines(beamStr)
  f.writelines(['!========================== Lattice ========================== \n'])
  f.writelines(latticeStr)
  f.close()


def _readraw(fname):
  try:
    fin = open(fname,'r')
    dataList  = fin.readlines()
    fin.close()
  except:
    print(( "  ERROR! Can't open file '" + fname + "'"))
    return False

  i=0
  while i < len(dataList):
    if dataList[i].lstrip()=='' or dataList[i].lstrip().startswith('!'):
      del dataList[i]
      i=i-1
    else:
      index = dataList[i].lstrip().find('!')
      if index==-1:
          dataList[i]=dataList[i].strip()+'\n'
      else:
          dataList[i]=dataList[i].lstrip()[:index].rstrip()+'\n'
    i=i+1
  dataList  = [line.split() for line in dataList ]
  for i in range(0,len(dataList)):
    for j in range(0,len(dataList[i])):
      if dataList[i][j] == '/':
        del dataList[i][j:]
    for j in range(0,len(dataList[i])):
      dataList[i][j] = _DtoE(dataList[i][j])
  return dataList

def _DtoE(word):
  if 'D' in word or 'd' in word: 
    try:
      temp = float(word.replace('D','E',1).replace('d','e',1))
      return str(temp)
    except:
      return word
  else:
    return word


def _str2beam(raw):
  """
  beam = _str2beam(beamStr)
  
  from (IMPACT format) string to beam class
  input 
      raw = (str) split string list of a IMPACT input
  output 
      beam = (dict) beam class
  """
  '''cpu cores'''
  i=0
  beam = data.beam({})
  print('  : mpi task info .......................',end='')
  beam.nCore_y = int(raw[i][0])
  beam.nCore_z = int(raw[i][1])
  '''simulation control parameters'''
  i+=1
  print('......done')
  print('  : simulation control parameters .......',end='')
  beam.dim             = int(raw[i][0])
  beam.n_particles     = int(raw[i][1])
  beam.integrator      = data.integrator[int(raw[i][2])]
  beam.error_study     = int(raw[i][3])==1
  beam.standard_output = data.standard_output[int(raw[i][4])]
  '''space charge field solver mesh info'''
  i+=1
  print('......done')
  print('  : space charge field solver, mesh info ',end='')
  mesh = data.dictClass({'mesh_x':int(raw[i][0]),
                    'mesh_y':int(raw[i][1]),
                    'mesh_z':int(raw[i][2]),
                    'fld_solver':data.fld_solver[int(raw[i][3])],
                    'boundary_x':float(raw[i][4]),
                    'boundary_y':float(raw[i][5]),
                    'boundary_z':float(raw[i][6]) } )
  beam.mesh = mesh
  '''Distribution'''
  i+=1
  print('......done')
  print('  : dist-type,restart,subcycle,#of state ',end='')
  distribution = data.dictClass( {'distribution_type':data.distribution_type[int(raw[i][0])]} )
  beam.restart = int(raw[i][1])==1
  beam.subcycle= int(raw[i][2])==1
  multi_charge = data.dictClass({'n_states':int(raw[i][3])})
  '''Multiple Charge State'''
  i+=1
  print('......done')
  print('  : Multiple Charge State info ..........',end='')
  multi_charge.n_particles = []
  for j in range(multi_charge.n_states):
    multi_charge.n_particles.append(int(raw[i][j]))
  i+=1
  multi_charge.current = []
  for j in range(multi_charge.n_states):
    multi_charge.current.append(float(raw[i][j]))
  i+=1
  multi_charge.q_m = []
  for j in range(multi_charge.n_states):
    multi_charge.q_m.append(float(raw[i][j]))
  beam.multi_charge = multi_charge
  '''Twiss'''
  i+=1
  print('......done')
  print('  : particle distribution info ..........',end='')
  if distribution.distribution_type in ['ReadFile_binary','IOTA_Waterbag','IOTA_Gauss']:
    if distribution.distribution_type == 'ReadFile_binary':
      try:
        distribution.file_id  = float(raw[i][0])
      except:
        print('file_id need to be provided in beam.distribution')
    else:
      distribution.NL_t  = float(raw[i][0])
      distribution.NL_c  = float(raw[i][1])
      distribution.betx  = float(raw[i][2])
      distribution.alfx  = -float(raw[i][3])/2.0
      distribution.emitx = float(raw[i][4])
      if distribution.distribution_type == 'IOTA_Gauss':
        distribution.CL    = float(raw[i][5])
    i+=1
    distribution.offsetx  = float(raw[i][0])
    distribution.offsetpx = float(raw[i][1])
    distribution.offsety  = float(raw[i][2])
    distribution.offsetpy = float(raw[i][3])
    i+=1
  else:
    distribution.sigmax  = float(raw[i][0])
    distribution.lambdax = float(raw[i][1])
    distribution.mux     = float(raw[i][2])
    distribution.scalex  = float(raw[i][3])
    distribution.scalepx = float(raw[i][4])
    distribution.offsetx = float(raw[i][5])
    distribution.offsetpx= float(raw[i][6])
    i+=1
    distribution.sigmay  = float(raw[i][0])
    distribution.lambday = float(raw[i][1])
    distribution.muy     = float(raw[i][2])
    distribution.scaley  = float(raw[i][3])
    distribution.scalepy = float(raw[i][4])
    distribution.offsety = float(raw[i][5])
    distribution.offsetpy= float(raw[i][6])
    i+=1
  distribution.sigmaz  = float(raw[i][0])
  distribution.lambdaz = float(raw[i][1])
  distribution.muz     = float(raw[i][2])
  distribution.scalez  = float(raw[i][3])
  distribution.scalepz = float(raw[i][4])
  distribution.offsetz = float(raw[i][5])
  distribution.offsetpz= float(raw[i][6])
  distribution.mode = 'impactdist'
  """ Reference Orbit """
  print('......done')
  print('  : beam reference orbit info ...........',end='')
  i+=1
  beam.current        = float(raw[i][0])
  beam.kinetic_energy = float(raw[i][1])
  beam.mass           = float(raw[i][2])
  beam.charge         = float(raw[i][3])
  beam.frequency      = float(raw[i][4])
  beam.phase          = float(raw[i][5])
  """ convert impact dist parameters to twiss parameters"""
  gamma = 1.0 + beam.kinetic_energy/beam.mass
  print('......done')
  print('  : converting impact dist to twiss param',end='')
  beam.distribution = distribution 
  beam.impactdist2twiss()
  print('......done')
  return beam

def _beam2str(beam):
  """
  beamStr = _beam2str(beam)
  
  from beam class class to (IMPACT format) string
  input 
      beam = beam class
  output 
      beamStr = split list of strings
  """
  '''cpu cores'''
  temp = [beam.nCore_y,beam.nCore_z]
  beamStr = [temp]
  
  '''simulation control parameters'''
  temp = [beam.dim,beam.n_particles,
          data.integrator.find_key(beam.integrator),
          int(beam.error_study),
          data.standard_output.find_key(beam.standard_output)]
  beamStr.append(temp)
  
  '''space charge field solver mesh info'''
  temp = [beam.mesh.mesh_x,beam.mesh.mesh_y,beam.mesh.mesh_z,
          data.fld_solver.find_key(beam.mesh.fld_solver),
          beam.mesh.boundary_x,beam.mesh.boundary_y,beam.mesh.boundary_z]
  beamStr.append(temp)
  
  '''Distribution'''
  temp = [data.distribution_type.find_key(beam.distribution.distribution_type),
          int(beam.restart),int(beam.subcycle),beam.multi_charge.n_states]
  beamStr.append(temp)
  
  '''Multiple Charge State'''
  beamStr.append(copy(beam.multi_charge.n_particles))
  beamStr.append(copy(beam.multi_charge.current))
  beamStr.append(copy(beam.multi_charge.q_m))
  
  '''Twiss'''
  gamma = 1.0 + beam.kinetic_energy/beam.mass
  #distribution = _twiss2impactdist(gamma,beam.frequency,beam.mass,beam.distribution)
  beam.twiss2impactdist()
  distribution = beam.distribution
  if distribution.distribution_type in ['ReadFile_binary','IOTA_Waterbag','IOTA_Gauss']:
    if distribution.distribution_type == 'ReadFile_binary':
      temp = [distribution.file_id,0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    else:
      temp = [distribution.NL_t,
              distribution.NL_c,
              distribution.betx,
              -2.0*distribution.alfx,
              distribution.emitx]
      if distribution.distribution_type == 'IOTA_Gauss':
        temp.append(distribution.CL)
      else:
        temp.append(0.0)
      temp.append(0.0)
    beamStr.append(temp)
    temp = [distribution.offsetx,
            distribution.offsetpx,
            distribution.offsety,
            distribution.offsetpy,
            0.0,0.0,0.0]
    beamStr.append(temp)
    temp = [distribution.sigmaz,
            distribution.lambdaz,
            distribution.muz,
            distribution.scalez,
            distribution.scalepz,
            distribution.offsetz,
            distribution.offsetpz]
    beamStr.append(temp)
  else:
    temp = [distribution.sigmax,
            distribution.lambdax,
            distribution.mux,
            distribution.scalex,
            distribution.scalepx,
            distribution.offsetx,
            distribution.offsetpx]
    beamStr.append(temp)
    temp = [distribution.sigmay,
            distribution.lambday,
            distribution.muy,
            distribution.scaley,
            distribution.scalepy,
            distribution.offsety,
            distribution.offsetpy]
    beamStr.append(temp)
    temp = [distribution.sigmaz,
            distribution.lambdaz,
            distribution.muz,
            distribution.scalez,
            distribution.scalepz,
            distribution.offsetz,
            distribution.offsetpz]
    beamStr.append(temp)
    
  """ Reference Orbit """
  temp = [beam.current,
          beam.kinetic_energy,
          beam.mass,
          beam.charge,
          beam.frequency,
          beam.phase]
  beamStr.append(temp)
  for i in range(len(beamStr)):
    for j in range(len(beamStr[i])):
      beamStr[i][j] = str(beamStr[i][j])
  return beamStr


def _str2elem(elemStr): 
  """
  elemDict = str2elem(elemStr)
  
  from (IMPACT format) string to element  
  input 
      elemStr = (str) string of a IMPACT lattice line
  output 
      elemDict = (dict) element dictionary 
  """
  elemID=intStr(elemStr[3])
  try:
    data.elem_type[elemID]
  except:
    print()
    print(data.bcolors.FAIL+'input element type code',elemID,
          ' is not listed in data.elem_type'+data.bcolors.END)
  if data.elem_type[elemID] == 'drift':
    elemDict = {'length'     : float(elemStr[0]),
                'n_sckick'   : intStr(elemStr[1]),
                'n_map'      : intStr(elemStr[2]),
                'pipe_radius': float(elemStr[4])
                }
  elif data.elem_type[elemID] == 'quad':
    elemDict = {'length'       : float(elemStr[0]),
                'n_sckick'     : intStr(elemStr[1]),
                'n_map'        : intStr(elemStr[2]),
                'B1'           : float(elemStr[4]),
                'file_id'      : intStr(elemStr[5]),
                'pipe_radius'  : float(elemStr[6])}
    if len(elemStr)>=8:
                 elemDict['misalign_x']=float(elemStr[7])
    if len(elemStr)>=9:
                 elemDict['misalign_y']=float(elemStr[8])
    if len(elemStr)>=10:
                 elemDict['rotation_x']=float(elemStr[9])
    if len(elemStr)>=11:
                 elemDict['rotation_y']=float(elemStr[10])
    if len(elemStr)>=12:
                 elemDict['rotation_z']=float(elemStr[11])
                 
  elif data.elem_type[elemID] == 'const_focusing':
    elemDict = { 'length'  : float(elemStr[0]),
                 'n_sckick': intStr(elemStr[1]),
                 'n_map'   : intStr(elemStr[2]),
                 'kx2'     : float(elemStr[4]),
                 'ky2'     : float(elemStr[5]),
                 'kz2'     : float(elemStr[6]),
                 'pipe_radius': float(elemStr[7])}
                 
  elif data.elem_type[elemID] == 'solenoid':
    elemDict = { 'length'  : float(elemStr[0]),
                 'n_sckick': intStr(elemStr[1]),
                 'n_map'   : intStr(elemStr[2]),
                 'Bz'      : float(elemStr[4]),
                 'file_id' : intStr(elemStr[5]),
                 'pipe_radius': float(elemStr[7])}
    if len(elemStr)>=9:
                 elemDict['misalign_x']=float(elemStr[8])
    if len(elemStr)>=10:
                 elemDict['misalign_y']=float(elemStr[9])
    if len(elemStr)>=11:
                 elemDict['rotation_x']=float(elemStr[10])
    if len(elemStr)>=12:
                 elemDict['rotation_y']=float(elemStr[11])
    if len(elemStr)>=13:
                 elemDict['rotation_z']=float(elemStr[12])
                 
  elif data.elem_type[elemID] == 'dipole':
    elemDict = { 'length'       : float(elemStr[0]),
                 'n_sckick'     : intStr(elemStr[1]), 
                 'n_map'        : intStr(elemStr[2]), 
                 'bending_angle': float(elemStr[4]), 
                 'k1'           : float(elemStr[5]), 
                 'file_id'      : int(float(elemStr[6])),
                 'pipe_radius'  : float(elemStr[7])
                }
    if len(elemStr)>=9:
                 elemDict['entrance_angle']=float(elemStr[8])
    if len(elemStr)>=10:
                 elemDict['exit_angle']=float(elemStr[9])
    if len(elemStr)>=11:
                 elemDict['entrance_curvature']=float(elemStr[10])
    if len(elemStr)>=12:
                 elemDict['exit_curvature']=float(elemStr[11])
    if len(elemStr)>=13:
                 elemDict['fringe_field_integration']=float(elemStr[12])

  elif data.elem_type[elemID] == 'multipole_thin':
    elemDict = {'KL_dipole': float(elemStr[5]),
                'KL_quad'  : float(elemStr[6]),
                'KL_sext'  : float(elemStr[7])}
    if len(elemStr)>=9:
                 elemDict['KL_oct']=float(elemStr[8])
    if len(elemStr)>=10:
                 elemDict['KL_deca']=float(elemStr[9])
    if len(elemStr)>=11:
                 elemDict['KL_dodeca']=float(elemStr[10])
 
  elif data.elem_type[elemID] == 'linear_matrix_map':
    elemDict = {
                 'nonlinear_insert_length' : float(elemStr[5]), 
                 'nonlinear_insert_tuneAdvance': float(elemStr[6]), 
                 'tune_advance' : float(elemStr[7]),
                }

  elif data.elem_type[elemID] in ['nonlinear_insert','nonlinear_insert_smooth_focusing']:
    elemDict = { 'length'            : float(elemStr[0]),
                 'n_sckick'          : intStr(elemStr[1]), 
                 'n_map'             : intStr(elemStr[2]), 
                 'strength_t'        : float(elemStr[4]), 
                 'transverse_scale_c': float(elemStr[5]), 
                 'pipe_radius'       : float(elemStr[7])
                }
    if data.elem_type[elemID] == 'nonlinear_insert':
       elemDict['tune_advance'] = float(elemStr[6])
    else:
      elemDict['betx'] = float(elemStr[6])
      
  elif data.elem_type[elemID] == 'DTL':
    elemDict= { 'length' : float(elemStr[0]),
                'n_sckick': intStr(elemStr[1]), 
                'n_map': intStr(elemStr[2]), 
                'field_scaling': float(elemStr[4]), 
                'frequency': float(elemStr[5]), 
                'phase': float(elemStr[6]), 
                'file_id': int(float(elemStr[7])),
                'pipe_radius': float(elemStr[8]),
                'quad1_length': float(elemStr[9]),
                'quad1_B1': float(elemStr[10]),
                'quad2_length': float(elemStr[11]),
                'quad2_B1': float(elemStr[12])
                }
    if len(elemStr)>=14:
                 elemDict['misalign_x']=float(elemStr[13])
    if len(elemStr)>=15:
                 elemDict['misalign_y']=float(elemStr[14])

  elif data.elem_type[elemID] == 'loop':
    elemDict = {'turns' : int(float(elemStr[5]))}
    
  elif data.elem_type[elemID] in ['CCDTL','CCL','SCRF','solenoidRF','EMfld']:
    elemDict= { 'length' : float(elemStr[0]),
                'n_sckick': intStr(elemStr[1]), 
                'n_map': intStr(elemStr[2]), 
                'field_scaling': float(elemStr[4]),
                'frequency': float(elemStr[5]),
                'phase': float(elemStr[6]),
                'file_id': int(float(elemStr[7])),
                'pipe_radius': float(elemStr[8])
                }
    if len(elemStr)>=10:
                 elemDict['misalign_x']=float(elemStr[9])
    if len(elemStr)>=11:
                 elemDict['misalign_y']=float(elemStr[10])
    if len(elemStr)>=12:
                 elemDict['rotation_x']=float(elemStr[11])
    if len(elemStr)>=13:
                 elemDict['rotation_y']=float(elemStr[12])
    if len(elemStr)>=14:
                 elemDict['rotation_z']=float(elemStr[13])
    if data.elem_type[elemID] == 'solenoidRF':
                 elemDict['Bz']=float(elemStr[14])

  elif data.elem_type[elemID] == 'centroid_shift':
    elemDict= {'x' : float(elemStr[5]),
               'px': float(elemStr[6]),
               'y' : float(elemStr[7]),
               'py': float(elemStr[8]),
               'z' : float(elemStr[9]),
               'pz': float(elemStr[10])}
    
  
  elif data.elem_type[elemID] == '-8':
    elemDict= {'file_id': intStr(elemStr[2]),
               'value'  : intStr(elemStr[4])}

  elif data.elem_type[elemID] == 'write_raw_ptcl':
    elemDict= {'file_id'  : intStr(elemStr[2]),
               'format_id': intStr(elemStr[4]),
               'turn'     : intStr(elemStr[5])}
    if len(elemStr)>=7:
      elemDict['sample_period']=intStr(elemStr[6])

  elif data.elem_type[elemID] == 'pipe_override':
    elemDict= {'pipe_shape': data.pipe_shape[intStr(elemStr[4])],
               'xmax'   : intStr(elemStr[5]),
               'ymax'   : intStr(elemStr[6])}

  elif data.elem_type[elemID] == 'TBT':
    elemDict= {'file_id'  : intStr(elemStr[2]),
               'pID_begin': intStr(elemStr[4]),
               'pID_end'  : intStr(elemStr[5])}
    
  elif data.elem_type[elemID] in ['TBT_integral','TBT_integral_onMomentum']:
    elemDict= {'file_id'           : intStr(elemStr[2]),
               'betx'              : float(elemStr[4]),
               'alfx'              : float(elemStr[5]),
               'strength_t'        : float(elemStr[6]), 
               'transverse_scale_c': float(elemStr[7]),
               'pID_begin'         : intStr(elemStr[8]),
               'pID_end'           : intStr(elemStr[9])}
  else :
    elemDict= {}
  elemDict['type']   = data.elem_type[elemID]
  elemDict = data.dictClass(elemDict)
  return elemDict
  
def _elem2str(elemDict): 
  """
  elemStr = str2elem(elemDict)
  
  from (IMPACT format) element to string
  input
    elemDict = (dict) element dictionary 
  output 
    elemStr = (str) string list of a IMPACT lattice line
  """
  try:
    elemStr = [elemDict.length, 0, 0, data.elem_type.find_key(elemDict.type)]
  except:
    elemStr = [0.0, 0, 0, data.elem_type.find_key(elemDict.type)]
  if elemDict.type == 'drift':
    elemStr[1]=elemDict.n_sckick
    elemStr[2]=elemDict.n_map
    elemStr.append(elemDict.pipe_radius)

  elif elemDict.type == 'quad':
    elemStr[1]=elemDict.n_sckick
    elemStr[2]=elemDict.n_map
    elemStr.append(elemDict.B1)
    elemStr.append(elemDict.file_id)
    elemStr.append(elemDict.pipe_radius)
    if 'misalign_x' in elemDict:
      elemStr.append(elemDict.misalign_x)
      if 'misalign_y' in elemDict:
        elemStr.append(elemDict.misalign_y)
        if 'rotation_x' in elemDict:
          elemStr.append(elemDict.rotation_x)
          if 'rotation_y' in elemDict:
            elemStr.append(elemDict.rotation_y)
            if 'rotation_z' in elemDict:
              elemStr.append(elemDict.rotation_z)
              
  elif elemDict.type == 'const_focusing':
    elemStr[1]=elemDict.n_sckick
    elemStr[2]=elemDict.n_map
    elemStr.append(elemDict.kx2)
    elemStr.append(elemDict.ky2)
    elemStr.append(elemDict.kz2)
    elemStr.append(elemDict.pipe_radius)

  elif elemDict.type == 'solenoid':
    elemStr[1]=elemDict.n_sckick
    elemStr[2]=elemDict.n_map
    elemStr.append(elemDict.Bz)
    elemStr.append(elemDict.file_id)
    elemStr.append(elemDict.pipe_radius)
    if 'misalign_x' in elemDict:
      elemStr.append(elemDict.misalign_x)
      if 'misalign_y' in elemDict:
        elemStr.append(elemDict.misalign_y)
        if 'rotation_x' in elemDict:
          elemStr.append(elemDict.rotation_x)
          if 'rotation_y' in elemDict:
            elemStr.append(elemDict.rotation_y)
            if 'rotation_z' in elemDict:
              elemStr.append(elemDict.rotation_z)
                 
  elif elemDict.type == 'dipole':
    elemStr[1]=elemDict.n_sckick
    elemStr[2]=elemDict.n_map
    elemStr.append(elemDict.bending_angle)
    elemStr.append(elemDict.k1)
    elemStr.append(elemDict.file_id)
    elemStr.append(elemDict.pipe_radius)
    if 'entrance_angle' in elemDict:
      elemStr.append(elemDict.entrance_angle)
      if 'exit_angle' in elemDict:
        elemStr.append(elemDict.exit_angle)
        if 'entrance_curvature' in elemDict:
          elemStr.append(elemDict.entrance_curvature)
          if 'exit_curvature' in elemDict:
            elemStr.append(elemDict.exit_curvature)
            if 'fringe_field_integration' in elemDict:
              elemStr.append(elemDict.fringe_field_integration)

  elif elemDict.type == 'multipole_thin':
    elemStr[0]=0.0
    elemStr.append(0.0)
    elemStr.append(elemDict.KL_dipole)
    elemStr.append(elemDict.KL_quad)
    elemStr.append(elemDict.KL_sext)
    if 'KL_oct' in elemDict:
      elemStr.append(elemDict.KL_oct)
      if 'KL_deca' in elemDict:
        elemStr.append(elemDict.KL_deca)
        if 'KL_dodeca' in elemDict:
          elemStr.append(elemDict.KL_dodeca)


  elif elemDict.type == 'linear_matrix_map':
    elemStr.append(0.0)
    elemStr.append(elemDict.nonlinear_insert_length)
    elemStr.append(elemDict.nonlinear_insert_tuneAdvance)
    elemStr.append(elemDict.tune_advance)

  elif elemDict.type in ['nonlinear_insert','nonlinear_insert_smooth_focusing']:
    elemStr[1]=elemDict.n_sckick
    elemStr[2]=elemDict.n_map
    elemStr.append(elemDict.strength_t)
    elemStr.append(elemDict.transverse_scale_c)
    if elemDict.type == 'nonlinear_insert':
      elemStr.append(elemDict.tune_advance)
    else:
      elemStr.append(elemDict.betx)
    elemStr.append(elemDict.pipe_radius)

  elif elemDict.type == 'DTL':
    elemStr[1]=elemDict.n_sckick
    elemStr[2]=elemDict.n_map
    elemStr.append(elemDict.field_scaling)
    elemStr.append(elemDict.frequency)
    elemStr.append(elemDict.phase)
    elemStr.append(elemDict.file_id)
    elemStr.append(elemDict.pipe_radius)
    elemStr.append(elemDict.quad1_length)
    elemStr.append(elemDict.quad1_B1)
    elemStr.append(elemDict.quad2_length)
    elemStr.append(elemDict.quad2_B1)
    if 'misalign_x' in elemDict:
      elemStr.append(elemDict.misalign_x)
      if 'misalign_y' in elemDict:
        elemStr.append(elemDict.misalign_y)

  elif elemDict.type == 'loop':
    elemStr.append(0.0)
    elemStr.append(elemDict.turns)
    
  elif elemDict.type in ['CCDTL','CCL','SCRF','solenoidRF','EMfld']:
    elemStr[1]=elemDict.n_sckick
    elemStr[2]=elemDict.n_map
    elemStr.append(elemDict.field_scaling)
    elemStr.append(elemDict.frequency)
    elemStr.append(elemDict.phase)
    elemStr.append(elemDict.file_id)
    elemStr.append(elemDict.pipe_radius)
    if 'misalign_x' in elemDict:
      elemStr.append(elemDict.misalign_x)
      if 'misalign_y' in elemDict:
        elemStr.append(elemDict.misalign_y)
        if 'rotation_x' in elemDict:
          elemStr.append(elemDict.rotation_x)
          if 'rotation_y' in elemDict:
            elemStr.append(elemDict.rotation_y)
            if 'rotation_z' in elemDict:
              elemStr.append(elemDict.rotation_z)
              if elemDict.type == 'solenoidRF':
                elemStr.append(elemDict.Bz)

  elif elemDict.type == 'centroid_shift':
    elemStr.append(1.0)
    elemStr.append(elemDict.x)
    elemStr.append(elemDict.px)
    elemStr.append(elemDict.y)
    elemStr.append(elemDict.py)
    elemStr.append(elemDict.z)
    elemStr.append(elemDict.pz)
    
  elif elemDict.type == '-8':
    elemStr[2]=elemDict.file_id
    elemStr.append(elemDict.value)
    
  elif elemDict.type == 'write_raw_ptcl':
    elemStr[2]=elemDict.file_id
    elemStr.append(elemDict.format_id)
    elemStr.append(elemDict.turn)
    if 'sample_period' in elemDict:
      elemStr.append(elemDict.sample_period)
    
  elif elemDict.type == 'TBT':
    elemStr[2]=elemDict.file_id
    elemStr.append(elemDict.pID_begin)
    elemStr.append(elemDict.pID_end)

  elif elemDict.type == 'pipe_override':
    elemStr.append(data.pipe_shape.find_key(elemDict.pipe_shape))
    elemStr.append(elemDict.xmax)
    elemStr.append(elemDict.ymax)
    
  elif elemDict.type in ['TBT_integral','TBT_integral_onMomentum']:
    elemStr[2]=elemDict.file_id
    elemStr.append(elemDict.betx)
    elemStr.append(elemDict.alfx)
    elemStr.append(elemDict.strength_t)
    elemStr.append(elemDict.transverse_scale_c)
    elemStr.append(elemDict.pID_begin)
    elemStr.append(elemDict.pID_end)
                
  for i in range(len(elemStr)):
    elemStr[i] = str(elemStr[i])

  return elemStr

  
def readReferenceOrbit(fileloc=''):
  file = open(fileloc+'fort.18','r')
  lines = file.readlines()
  file.close()
  f=data.dictClass({'s':[],'phase':[],'gamma':[],
               'beta':[],'max_radius':[]})
  for j in range(len(lines)) :
    lines[j] = lines[j].split()
    f.s.append(float(lines[j][0]))
    f.phase.append(float(lines[j][1]))
    f.gamma.append(float(lines[j][2]))
    f.beta.append(float(lines[j][4]))
    f.max_radius.append(float(lines[j][5]))
  return f

def readReferenceOrbitAt(sIndex,fileloc=''):
  file = open(fileloc+'fort.18','r')
  lines = file.readlines()
  lines = [float(lines[sIndex].split()[i]) for i in range(5)]
  file.close()
  f=data.dictClass({'s':lines[0],'phase':lines[1],
               'gamma':lines[2],'beta':lines[4],'max_radius':lines[5]})
  return f

  
def get_sIndex(s,fileloc=''):
  rf = readReferenceOrbit(fileloc)
  for i in range(1,len(rf.s)) :
    s1 =  rf.s[i]
    if s-s1 < 0 :
      break
  if s1-s < s-rf.s[i-1]:
    return i,rf.s[i]
  else:
    return i-1,rf.s[i-1]

    
def readRMS(direction, sample_period=1,fileLoc=''):
  """
  f = readRMS(direction,sample_period=1,fileLoc='')
  Read RMS beam info
  input 
    direction = (char) 'x', 'y' or 'z'
    sample_period = (int>0) number of lines to skip when reading output    
  """
  if direction == 'x':
    file = open(fileLoc+'fort.24','r')
    lines = file.readlines()
    file.close()
    f=data.dictClass({'s':[],'centroid_x':[],'rms_x':[],
                 'centroid_px':[],'rms_px':[],
                 'alfx':[],'emitx':[]})
    for j in range(0,len(lines),sample_period) :
      lines[j] = lines[j].split()
      f.s.append(float(lines[j][0]))
      f.centroid_x.append(float(lines[j][1]))
      f.rms_x.append(float(lines[j][2]))
      f.centroid_px.append(float(lines[j][3]))
      f.rms_px.append(float(lines[j][4]))
      f.alfx.append(float(lines[j][5]))
      f.emitx.append(float(lines[j][6]))
      
  elif direction == 'y':
    file = open(fileLoc+'fort.25','r')
    lines = file.readlines()
    file.close()
    f=data.dictClass({'s':[],'centroid_y':[],'rms_y':[],
                 'centroid_py':[],'rms_py':[],
                 'alfy':[],'emity':[]})
    for j in range(0,len(lines),sample_period) :
      lines[j] = lines[j].split()
      f.s.append(float(lines[j][0]))
      f.centroid_y.append(float(lines[j][1]))
      f.rms_y.append(float(lines[j][2]))
      f.centroid_py.append(float(lines[j][3]))
      f.rms_py.append(float(lines[j][4]))
      f.alfy.append(float(lines[j][5]))
      f.emity.append(float(lines[j][6]))
      
  elif direction == 's':
    file = open(fileLoc+'fort.26','r')
    lines = file.readlines()
    file.close()
    f=data.dictClass({'s':[],'centroid_z':[],'rms_z':[],
                 'centroid_pz':[],'rms_pz':[],
                 'alfz':[],'emitz':[]})
    for j in range(0,len(lines),sample_period) :
      lines[j] = lines[j].split()
      f.s.append(float(lines[j][0]))
      f.centroid_z.append(float(lines[j][1]))
      f.rms_z.append(float(lines[j][2]))
      f.centroid_pz.append(float(lines[j][3]))
      f.rms_pz.append(float(lines[j][4]))
      f.alfz.append(float(lines[j][5]))
      f.emitz.append(float(lines[j][6])*1.0e6) # degree-MeV
  for k in f.keys():
      f[k] = np.array(f[k])
  return f
  
def readRMS_at(sIndex,direction,fileLoc=''):
  """
  f = readRMS_at(sIndex,direction,sample_period=1,fileLoc=''):
  Read RMS beam size at location corresponds to sIndex
  input 
    direction = (char) 'x', 'y' or 'z'
    sample_period = (int>0) number of lines to skip when reading output    
  """
  if direction == 'x':
    file = open(fileLoc+'fort.24','r')
    lines = file.readlines()
    lines = lines[sIndex].split()
    file.close()
    f=data.dictClass({'s':float(lines[0]),
                 'centroid_x':float(lines[1]),
                 'rms_x':float(lines[2]),
                 'centroid_px':float(lines[3]),
                 'rms_px':float(lines[4]),
                 'alfx':float(lines[5]),
                 'emitx':float(lines[6])})
  elif direction == 'y':
    file = open(fileLoc+'fort.25','r')
    lines = file.readlines()
    lines = lines[sIndex].split()
    file.close()
    f=data.dictClass({'s':float(lines[0]),
                 'centroid_y':float(lines[1]),
                 'rms_y':float(lines[2]),
                 'centroid_py':float(lines[3]),
                 'rms_py':float(lines[4]),
                 'alfy':float(lines[5]),
                 'emity':float(lines[6])})
      
  elif direction == 'z':
    file = open(fileLoc+'fort.26','r')
    lines = file.readlines()
    lines = lines[sIndex].split()
    file.close()
    f=data.dictClass({'s':float(lines[0]),
                 'centroid_z':float(lines[1]),
                 'rms_z':float(lines[2]),
                 'centroid_pz':float(lines[3]), 
                 'rms_pz':float(lines[4]),
                 'alfz':float(lines[5]),
                 'emitz':float(lines[6]) # degree-MeV
                })
  return f


def readOptics(direction,sample_period=1,fileLoc=''):
  """
  f = readOptics(direction,sample_period=1,fileLoc='')
  Read Optics functions ( Optics ftn is calcualted using beam porfile)
  *** note that it is not optics of the lattice. 
      it is optics of the beam.
      Can be very different from lattice optics 
      when large mismatch or space charge is present. ****
  input 
      direction = 'x' or 'y' or 'z'
      sample_period  = sampling period from beam distribution output file of IMPACTz
      fileLoc = (string) path
  output 
      f = optics parameters at every (sample_period) integration step.
  """
  if direction == 'x':
    file = open(fileLoc+'fort.24','r')
    f=data.dictClass({'s':[],'betx':[],'alfx':[],'emitx':[],'phx':[]})
  elif direction == 'y':
    file = open(fileLoc+'fort.25','r')
    f=data.dictClass({'s':[],'bety':[],'alfy':[],'emity':[],'phy':[]})
  elif direction == 'z':
    file = open(fileLoc+'fort.26','r')
    f=data.dictClass({'s':[],'betz':[],'alfz':[],'emitz':[],'phz':[]})
  lines = file.readlines()
  file.close()
  for i in range(len(lines)):
    lines[i] = [ float(lines[i].split()[j]) for j in [0,2,4,5,6] ]
  ph = 0
  s0 = 0
  if direction == 'x':
    j=sample_period-1
    for i in range(len(lines)):
      s, sigmax, sigmap, alpha, emittance_norm = lines[i]
      beta = (1+alpha*alpha)**0.5 *sigmax/sigmap
      ph = ph + (s-s0)/beta
      s0 = s
      j+=1
      if j==sample_period:
        j=0
        f.s.append(s)
        f.betx.append(beta)
        f.alfx.append(alpha)
        f.emitx.append(emittance_norm)
        f.phx.append(ph)
  elif direction == 'y':
    j=sample_period-1
    for i in range(len(lines)):
      s, sigmax, sigmap, alpha, emittance_norm = lines[i]
      beta = (1+alpha*alpha)**0.5 *sigmax/sigmap
      ph = ph + (s-s0)/beta
      s0 = s
      j+=1
      if j==sample_period:
        j=0
        f.s.append(s)
        f.bety.append(beta)
        f.alfy.append(alpha)
        f.emity.append(emittance_norm)
        f.phy.append(ph)
  elif direction == 'z':
    j=sample_period-1
    for i in range(len(lines)):
      s, sigmax, sigmap, alpha, emittance_norm = lines[i]
      beta = (1+alpha*alpha)**0.5 *sigmax/sigmap*1.0e-6
      ph = ph + (s-s0)/beta
      s0 = s
      j+=1
      if j==sample_period:
        j=0
        f.s.append(s)
        f.betz.append(beta)
        f.alfz.append(alpha)
        f.emitz.append(emittance_norm*1.0e6)
        f.phz.append(ph)
  return f
  
def readOpticsAt(sIndex,direction,fileLoc=''):
  """
  f = readOpticsAt(sIndex,direction,fileLoc='')
  Read Optics functions ( Optics ftn is calcualted using beam porfile)
  *** note that it is not optics of the lattice. 
      it is optics of the beam.
      Can be very different from lattice optics 
      when large mismatch or space charge is present. ****
  input 
      direction = 'x' or 'y' or 'z'
      fileLoc = (string) path
  """
  f = readOptics(direction,fileLoc=fileLoc)
  for k,v in f.items():
    f[k] = v[sIndex]
  return f

    
def readLost(sample_period=1,fileLoc=''):
  file = open(fileLoc+'fort.32','r')
  lines = file.readlines()
  file.close()
  f=[]
  for i in range(0,len(lines),sample_period) :
    f.append( int(lines[i].split()[1]) )
  return f 
  
def readLostAt(sIndex, fileLoc=''):
  file = open(fileLoc+'fort.32','r')
  lines = file.readlines()
  file.close()
  return int(lines[sIndex].split()[1])

def readLost_detail(fileLoc=''):
  tmp = np.loadtxt(fileLoc+'lost_partcl.data',skiprows=1)
  f = data.dictClass()
  f['z']=tmp[:,0]
  f['x']=tmp[:,1]
  f['y']=tmp[:,2]
  intvec = np.vectorize(int)
  f['particle_id']=intvec(tmp[:,3])
  f['n_lost']=len(f.z)
  f['elem_type']=(data.elem_type[int(tmp[i,4])] for i in range(f.n_lost))
  f['n-th elem']=int(tmp[:,5])
  return f

#%%############################################################################
###############################################################################
###                      Particle data Manipulator                          ###
###############################################################################
############################################################################### 
def normalizeParticleData(data, ke, mass, freq):
  gamma = ke/mass+1.0
  beta = np.sqrt(1.0-1.0/(gamma*gamma))
  x_norm = 2*freq*3.141592653589793/299792458
  px_norm = gamma*beta
  f = copy(data)
  f[:,0] = data[:,0]*x_norm
  f[:,1] = data[:,1]*px_norm
  f[:,2] = data[:,2]*x_norm
  f[:,3] = data[:,3]*px_norm
  f[:,4] = np.pi/180*data[:,4]
  f[:,5] = -data[:,5]/mass
  return f
    
def unNormalizeParticleData(data, ke, mass, freq):
  gamma = ke/mass+1.0
  beta = np.sqrt(1.0-1.0/(gamma*gamma))
  x_norm = 2*freq*3.141592653589793/299792458
  px_norm = gamma*beta
  f = copy(data)
  f[:,0] = data[:,0]/x_norm
  f[:,1] = data[:,1]/px_norm
  f[:,2] = data[:,2]/x_norm
  f[:,3] = data[:,3]/px_norm
  f[:,4] = 180/np.pi*data[:,4]
  f[:,5] = -mass*data[:,5]
  return f

def readParticleData(fileID, ke, mass, freq, format_id=0,fileLoc=''):
    if isinstance(fileID, str):
      data=np.loadtxt(fileID,skiprows=1)
    elif format_id=0:
      data=np.loadtxt(fileLoc+'fort.'+str(fileID))
      return unNormalizeParticleData(data, ke, mass, freq)
    elif format_id=1:
      cPath = os.getcwd()
      os.chdir(cPath+fileLoc)
      if not os.path.isfile('fort.'+str(-fileID)):
        print('can not find <fort.'+str(-fileID)+'> file')
        return 'file error'
      npt = _read_pdata.read_phasespace_size(fileID)
      data= np.transpose(_read_pdata.read_phasespace(fileID,npt))
      os.chdir(cPath)
      return unNormalizeParticleData(data, ke, mass, freq)
    elif format_id=2:
      return
    
def readParticleDataSliced(nSlice, fileID, ke, mass, freq, zSliced=True, fileLoc=''):
    """
    pData = readParticleDataSliced(nSlice, fileID, ke, mass, freq, 
                                   zSliced=True, fileLoc='')
    
    slice particle data into energy or longtudinal bin
    input : 
        nSlice : (int) number of slices
        fileID : (int) field ID to read IMPACT particle output
        ke : (real) kinetic energy of the IMPACT particle output
        freq : (real) reference frequency used to define 
                      longitudinal coordinate
        zSliced : (bool) if False, energy bin sliced
        fileLoc : (string) path to the IMPACT particle output
    output
        pData : (numpy arr) slized particle data
    """
    if fileID>0:
      data=np.loadtxt(fileLoc+'fort.'+str(fileID))
    else:
      cPath = os.getcwd()
      os.chdir(cPath+fileLoc)
      npt = _read_pdata.read_phasespace_size(fileID)
      data= np.transpose(_read_pdata.read_phasespace(fileID,npt))
      os.chdir(cPath)
      
    datatmp=unNormalizeParticleData(data, ke, mass, freq)
    
    f=[]    
    if zSliced:
        z_min = min(datatmp[:,4])
        z_max = max(datatmp[:,4])
        dz = (z_max-z_min)/float(nSlice)
        for i in range(nSlice):
            temp = []
            for j in range(len(datatmp[:,4])):
                if z_min + i*dz < datatmp[j,4] < z_min + (i+1)*dz :
                    temp.append(data[j,:])
            f.append(np.array(temp))
        return f
        
    else:
        ke_min = min(datatmp[:,5])
        ke_max = max(datatmp[:,5])
        dke = (ke_max-ke_min)/float(nSlice)
        for i in range(nSlice):
            temp = []
            for j in range(len(datatmp[:,5])):
                if ke_min + i*dke < datatmp[j,5] < ke_min + (i+1)*dke :
                    temp.append(data[j,:])
            f.append(np.array(temp))  
        return f

#%%############################################################################
###############################################################################
###                           Lattice Manipulator                           ###
###############################################################################
############################################################################### 
def writeParticleData(data, ke, mass, freq, fileLoc='',fname='partcl.data'):
    tmpdata=normalizeParticleData(data, ke, mass, freq)
    np.savetxt(fname,tmpdata,header=str(len(data))+' 0. 0.',comments='')
