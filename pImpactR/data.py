clight = 299792458  # m/s
pi = 3.141592653589793
twopi = 2.0*pi

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

class beam(dictClass) :
  def __init__(self,dist_type='Waterbag'):
    self.nCore_y=1
    self.nCore_z=1
    self.dim=6
    self.n_particles=2048
    self.integrator='Linear'
    self.error_study=False
    self.standard_output='standard output'
    self.restart=False
    self.subcycle=False
    self.current=0.0
    self.kinetic_energy=2.5e6
    self.mass=938.272045e6
    self.charge=1.0
    self.frequency=30.0e6
    self.phase=0.0
    self.mesh = dictClass({
                'mesh_x':32,
                'mesh_y':32,
                'mesh_z':32,
                'fld_solver':'Trans:open,  Longi:open',
                'boundary_x':0.014,
                'boundary_y':0.014,
                'boundary_z':0.1
                })
    if dist_type in ['IOTA_Waterbag','IOTA_Gauss']:
      self.distribution = dictClass( {
                'distribution_type':dist_type,
                'NL_t' :0.0,
                'NL_c' :0.0,
                'betx' :0.0,
                'betpx':0.0,
                'emitx':0.0,
                'CL'   :0.0,
                'scalez'  :0.0,
                'scalepz' :0.0,
                'offsetz' :0.0,
                'offsetpz':0.0
                })
    else:
      self.distribution = dictClass({
                'distribution_type':dist_type,
                'betx' :0.0,
                'alfx' :0.0,
                'emitx':0.0,
                'bety' :0.0,
                'alfy' :0.0,
                'emity':0.0,
                'betz' :0.0,
                'alfz' :0.0,
                'emitz':0.0,
                'scalex'  :0.0,
                'scalepx' :0.0,
                'offsetx' :0.0,
                'offsetpx':0.0,
                'scaley'  :0.0,
                'scalepy' :0.0,
                'offsety' :0.0,
                'offsetpy':0.0,
                'scalez'  :0.0,
                'scalepz' :0.0,
                'offsetz' :0.0,
                'offsetpz':0.0
                })
    self.distribution.mode = 'twiss'
    multi_charge = dictClass({
                'n_states'   :1,
                'n_particles':[2048],
                'current'    :[0.0],
                'q_m'        :[1.0657889726792521e-09]
                 })

  def twiss2impactdist(self):
    """
    convert twiss parameters to Impact distparam
    """
    if self.distribution.mode == 'impactdist' :
      return
    gamma = 1.0+self.kinetic_energy/self.mass
    freq  = self.frequency
    mass  = self.mass
    x_norm = clight/(twopi*freq)
    bg = (gamma**2-1.0)**0.5
    px_norm = 1.0/bg
    z_norm = 180.0/pi
    pz_norm= mass/1.0e6
    twiss = self.distribution
    param = dictClass({'distribution_type':twiss.distribution_type})
    param.mode = 'impactdist'
    
    if param.distribution_type in ['IOTA_Waterbag','IOTA_Gauss']:
      if all (k in twiss.keys() for k in ('NL_t','NL_c')):
        param.NL_t = twiss.NL_t
        param.NL_c = twiss.NL_c
      else:
        raise KeyError('NL_t and NL_c must be present in beam.distribution for '+ param.distribution_type +' dist_type')
      param.betx  = twiss.betx  
      if 'betPx' in twiss.keys():
        param.betPx = twiss.betPx
      else:
        param.betPx = -2.0*twiss.alfx
      param.emitx = twiss.emitx
      if param.distribution_type == 'IOTA_Gauss':
        if 'CL' in twiss.keys():
          param.CL  = twiss.CL
        else:
          param.CL = 3.0
    else:
      betx  = twiss.betx
      emitx = twiss.emitx
      alfx  = twiss.alfx
      if 0.0 in [betx,emitx]:
        param.sigmax = 0.0
        param.lambdax = 0.0
      else:
        param.sigmax = (betx*(emitx/bg)/(1.0+alfx*alfx))**0.5/x_norm
        param.lambdax = ((emitx/bg)/betx)**0.5/px_norm
      param.mux = alfx / (1.0+alfx*alfx)**0.5 
      param.scalex  = twiss.scalex  
      param.scalepx = twiss.scalepx 
      param.offsetx = twiss.offsetx/x_norm
      param.offsetpx= twiss.offsetpx/px_norm
      
      bety  = twiss.bety
      emity = twiss.emity
      alfy  = twiss.alfy
      if 0.0 in [bety,emity]:
        param.sigmay = 0.0
        param.lambday = 0.0
      else:
        param.sigmay = (bety*(emity/bg)/(1.0+alfy*alfy))**0.5/x_norm
        param.lambday = ((emity/bg)/bety)**0.5/px_norm
      param.muy = alfy / (1.0+alfy*alfy)**0.5 
      param.scaley  = twiss.scaley  
      param.scalepy = twiss.scalepy 
      param.offsety = twiss.offsety/x_norm
      param.offsetpy= twiss.offsetpy/px_norm
      
    betz  = twiss.betz
    emitz = twiss.emitz 
    alfz  = twiss.alfz
    if 0.0 in [betz,emitz]:
      param.sigmaz = 0.0
      param.lambdaz = 0.0
    else:
      param.sigmaz = ( betz*emitz/(1.0+alfz*alfz) )**0.5/z_norm
      param.lambdaz = (emitz/betz)**0.5/pz_norm
    param.muz = alfz / (1.0+alfz*alfz)**0.5 
    param.scalez  = twiss.scalez  
    param.scalepz = twiss.scalepz 
    param.offsetz = twiss.offsetz/z_norm
    param.offsetpz= twiss.offsetpz/pz_norm
      
    self.distribution = param

  def impactdist2twiss(self):
    """
    convert Impact distribution parameters to twiss parameters
    """
    if self.distribution.mode == 'twiss' :
      return
    gamma = 1.0+self.kinetic_energy/self.mass
    freq  = self.frequency
    mass  = self.mass
    x_norm = clight/(twopi*freq)
    bg = (gamma**2-1.0)**0.5
    px_norm = 1.0/bg
    z_norm = 180.0/pi
    pz_norm= mass/1e6
    
    param = self.distribution
    twiss = dictClass({'distribution_type':param.distribution_type})
    twiss.mode = 'twiss'
    # z
    z00 = param.sigmaz
    z11 = param.lambdaz
    z01 = param.muz
    alf = z01/(1.0-z01**2)**0.5
    if z00 == 0.0:
      if z11 != 0.0:
        z00 = 1.0e-14
    else:
      if z11 == 0.0:
        z11 = 1.0e-14
    emit= z00*z_norm*z11*pz_norm*(1.0+alf**2)**0.5
    if z11 == 0.0:
      beta = 0.0
    else:
      beta= emit/(z11*pz_norm)**2
    twiss.betz = beta # degree/MeV 
    twiss.alfz = alf
    twiss.emitz= emit  # degree-MeV 
    twiss.scalez = param.scalez
    twiss.scalepz= param.scalepz
    twiss.offsetz = param.offsetz*z_norm
    twiss.offsetpz= param.offsetpz*pz_norm 
      
    # x,y
    if param.distribution_type in ['IOTA_Waterbag','IOTA_Gauss']:
      twiss.NL_t  = param.NL_t  
      twiss.NL_c  = param.NL_c  
      twiss.betx  = param.betx  
      twiss.betPx = param.betPx 
      twiss.emitx = param.emitx
      if param.distribution_type == 'IOTA_Gauss':
        twiss.CL  = param.CL

    else:
      x00 = param.sigmax
      x11 = param.lambdax
      x01 = param.mux
      alf = x01/(1.0-x01**2)**0.5
      if x00 == 0.0:
        if x11 != 0.0:
          x00 = 1.0e-15
      else:
        if x11 == 0.0:
          x11 = 1.0e-15  
      emit= x00*x_norm*x11*(1.0+alf**2)**0.5
      if x11 == 0.0:
        beta = 0.0
      else:
        beta= emit*bg/x11**2
      twiss.betx = beta
      twiss.alfx = alf
      twiss.emitx= emit
      twiss.scalex = param.scalex
      twiss.scalepx= param.scalepx
      twiss.offsetx = param.offsetx*x_norm
      twiss.offsetpx= param.offsetpx*px_norm
      
      y00 = param.sigmay
      y11 = param.lambday
      y01 = param.muy
      alf = y01/(1.0-y01**2)**0.5
      if y00 == 0.0:
        if y11 != 0.0:
          y00 = 1.0e-15
      else:
        if y11 == 0.0:
          y11 = 1.0e-15  
      emit= y00*x_norm*y11*(1.0+alf**2)**0.5
      if y11 == 0.0:
        beta = 0.0
      else:
        beta= emit*bg/y11**2
      twiss.bety = beta
      twiss.alfy = alf
      twiss.emity= emit
      twiss.scaley = param.scaley
      twiss.scalepy= param.scalepy
      twiss.offsety = param.offsety*x_norm
      twiss.offsetpy= param.offsetpy*px_norm
    
    self.distribution = twiss


      
mass = dictClass({'electron'   :510998.9461,
             'proton'     :938272081.3})

distribution_type = dictClass({
                    1 :'Uniform'   ,
                    2 :'Gauss'     ,
                    3 :'Waterbag'  ,
                    4 :'SemiGauss' ,
                    5 :'KV'        ,
                    23:'ReadFile'  ,
                    81:'IOTA_Waterbag'   ,
                    82:'IOTA_Gauss',
                    16:'Multi-charge-state Waterbag',
                    17:'Multi-charge-state Gaussian'})

fld_solver = dictClass( {1:'Trans:open,  Longi:open'  ,
                    2:'Trans:open,  Longi:period',
                    3:'Trans:Round, Longi:open'  ,
                    4:'Trans:Round, Longi:perod' ,
                    5:'Trans:Rect,  Longi:open'  ,
                    6:'Trans:Rect,  Longi:perod' ,
                    7:'Symplectic_Spectral_2D'   ,
                    8:'PIC_2D'                   })

standard_output = dictClass({1:'standard output',
                        2:'90,95,99 emittance output'})
                    
integrator = dictClass({1:'Linear'   ,
                   2:'NonLinear'})

elem_type = dictClass({0  :'drift'         ,
                  1  :'quad'          ,
                  2  :'const_focusing',
                  3  :'solenoid'      ,
                  4  :'dipole'        ,
                  5  :'multipole_thin',
                  6  :'nonlinear_insert',
                  101:'DTL'           ,
                  102:'CCDTL'         ,
                  103:'CCL'           ,
                  104:'SCRF'          ,
                  105:'solenoidRF'    ,
                  106:'TWS'           ,
                  110:'EMfld'         ,
                  -1 :'centroid_to_0' ,
                  -2 :'write_raw_ptcl',
                  -7 :'save4restart'  ,
                  -8 :'-8'            ,
                  -16:'loop_through_lattice',
                  -21:'centroid_shift',
                  -46:'linear_matrix_map',
                  -88:'TBToutput_single_particle',  # turn-by-turn
                  -89:'TBToutput_multi_particles',  # turn-by-turn
                  -99:'halt'           })

unit = dictClass({'length'       :'m',
             'n_sckick'     :'1',
             'field_scaling':'1.0',
             'frequency'    :'Hz',
             'mass'         :'eV',
             'charge'       :'e',
             'current'      :'A',
             'phase'        :'rad',
             'B1'           :'T/m',
             'pipe_radius'  :'m',
             'n_map'        :'1',
             'bending_angle':'rad',
             'misalign_x'   :'m',
             'misalign_y'   :'m',
             'rotation_x'   :'rad',
             'rotation_y'   :'rad',
             'rotation_z'   :'rad',
             'boundary_x'   :'m',
             'boundary_y'   :'m',
             'boundary_z'   :'m',
             'Bz'           :'T',
             's'            :'m',
             'x'            :'m',
             'px'           :'1',
             'y'            :'m',
             'py'           :'1',
             'z'            :'degree',
             'pz'           :'MeV',
             'betx'         :'m',
             'bety'         :'m',
             'betz'         :'degree/MeV',
             'phx'          :'rad',
             'phy'          :'rad',
             'phz'          :'rad-m-MeV/degree',
             'emitx'        :'m-rad',
             'emity'        :'m-rad',
             'emitz'        :'degree-MeV',
             'scalex'       :'1.0',
             'scaley'       :'1.0',
             'scalez'       :'1.0',
             'scalepx'      :'1.0',
             'scalepy'      :'1.0',
             'scalepz'      :'1.0',
             'offsetx'      :'m',
             'offsety'      :'m',
             'offsetz'      :'degree',
             'offsetpx'     :'rad',
             'offsetpy'     :'rad',
             'offsetpz'     :'MeV',
             'max_radius'   :'m',
             'centroid_x'   :'m',
             'centroid_y'   :'m',
             'centroid_z'   :'degree',
             'centroid_px'  :'rad',
             'centroid_py'  :'rad',
             'centroid_pz'  :'MeV',
             'kinetic_energy'    :'eV',
             'entrance_angle'    :'rad',
             'exit_angle'        :'rad',
             'entrance_curvature':'rad',
             'exit_curvature'    :'rad',
             'fringe_field_integration' :'1'
             })

class bcolors:
  HEADER = '\033[95m'
  BLUE = '\033[94m'
  GREEN = '\033[92m'
  WARN = '\033[93m'
  FAIL = '\033[91m'
  END = '\033[0m'
  BOLD = '\033[1m'
  UNDERLINE = '\033[4m'