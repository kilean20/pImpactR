import numpy as np

try:
  from pyNaff import pynaff as naff
except:
  print('util.naff is not available')
#from scipy import optimize
#%%============================================================================
#                                 beam physics                              
#==============================================================================
Mp=938.2720813E6
Me=0.510998910E6
cLight=299792458.0

def RMSemittance(X,Px):
    x  = X-np.mean(X)
    px = Px-np.mean(Px)
    return np.sqrt(np.var(x)*np.var(px)-np.mean(x*px)**2)

def gamma(ke,mass):
    return ke/mass+1.0
    
def beta(ke,mass):
    rel_gam=ke/mass+1.0
    return np.sqrt(1.0-1.0/rel_gam/rel_gam)
    
def rigidity(ke,mass):
    rel_gam=ke/mass+1.0
    rel_bet=np.sqrt(1.0-1.0/rel_gam/rel_gam)
    return cLight/(mass*rel_bet*rel_gam)

def bendingRadius(B,ke,mass):
    return 1.0/B/rigidity(ke,mass)
    
def B1toK(b1,ke,mass):
    rel_gam=ke/mass+1.0
    rel_bet=np.sqrt(1.0-1.0/rel_gam/rel_gam)
    rigidity=cLight/(mass*rel_bet*rel_gam)
    return b1*rigidity  

def KtoB1(K,ke,mass):
    rel_gam=ke/mass+1.0
    rel_bet=np.sqrt(1.0-1.0/rel_gam/rel_gam)
    rigidity=cLight/(mass*rel_bet*rel_gam)
    return K/rigidity    
    
def RFphaseAdvancePerMeter(freq,ke,mass):
    rel_gam=ke/mass+1.0
    rel_bet=np.sqrt(1.0-1.0/rel_gam/rel_gam)
    return freq*2.0*np.pi/cLight/rel_bet

#%%============================================================================
#                                 matrix map                               
#==============================================================================
def Mdrift(L):
    return np.matrix([[1.0,  L ],
                      [0.0, 1.0]] )

def Mbend(L,rho):
    return np.matrix( [[np.cos(L/rho),      rho*np.sin(L/rho)],
                       [-np.sin(L/rho)/rho, np.cos(L/rho)    ]] )

def Mquad(L,K):
    if abs(K)<1E-10:
        drift = Mdrift(0.5*L)
        kick = np.matrix( [[1.0,    0.0],
                           [K*L, 1.0 ]] )
        return drift*kick*drift
    elif K>0:
        k=np.sqrt(K)
        return np.matrix( [[np.cos(k*L),    np.sin(k*L)/k],
                           [-k*np.sin(k*L), np.cos(k*L)  ]] )
    else:
        k=np.sqrt(-K)
        return np.matrix([[np.cosh(k*L),    np.sinh(k*L)/k],
                          [k*np.sinh(k*L),  np.cosh(k*L)  ]] )     

#%%============================================================================
#                                  optics                              
#==============================================================================
def sigma(bet,alf):
    return np.matrix( [[bet,  -alf],
                       [-alf, (1.0+alf*alf)/bet]] )

def sigma_next(Sigma,M):
    return M*Sigma*np.transpose(M)
    
def dispersion_next(Eta, M, Lbend=0.0, rho=0.0):
    if rho==0:
        return M*Eta
    else:
        return M*Eta + [[ rho*( 1.0- np.cos(Lbend/rho) ) ],
                        [ np.sin(Lbend/rho)       ]]    
    
    
#%%============================================================================
#                                  achromat                         
#==============================================================================
def update_lattice(lattice, ke=306.195E6, mass=Mp, flag_updateB1=False):
    for i in range(len(lattice)):
        if lattice[i]['type']=='bend':
            lattice[i]['rho'] = lattice[i]['length']/lattice[i]['angle']
            lattice[i]['Mx'] = Mbend(lattice[i]['length'], lattice[i]['rho'])
            lattice[i]['My'] = Mdrift(lattice[i]['length'])
        elif lattice[i]['type']=='quad':
            if not 'K1' in lattice[i]:
                lattice[i]['K1'] = B1toK(lattice[i]['B1'], ke, mass)
            lattice[i]['Mx'] = Mquad(lattice[i]['length'], lattice[i]['K1'])
            lattice[i]['My'] = Mquad(lattice[i]['length'], -lattice[i]['K1'])
            if flag_updateB1:
                lattice[i]['B1'] = KtoB1(lattice[i]['K1'], ke, mass)
        elif lattice[i]['type']=='drift':
            lattice[i]['Mx'] = Mdrift(lattice[i]['length'])
            lattice[i]['My'] = lattice[i]['Mx']
        else:
            lattice[i]['Mx'] = Mdrift(0.0)
            lattice[i]['My'] = lattice[i]['Mx']

def achromat_condition(lattice):
    M=np.matrix([ [1.,0.],[0.,1.] ] )
    Eta = np.matrix([ [0.],[0.] ] )
    for i in range(len(lattice)):
        if lattice[i]['type']=='bend':
            rho = lattice[i]['rho'] 
            L = lattice[i]['length']
            Eta[0,0] = Eta[0,0] - M[1,1]*rho*(1.0-np.cos(L/rho)) - M[0,1]*np.sin(L/rho)
            Eta[1,0] = Eta[1,0] + M[1,0]*rho*(1.0-np.cos(L/rho)) + M[0,0]*np.sin(L/rho)
        M=lattice[i]['Mx']*M
    Eta = M*Eta
        
    return Eta[0,0]**2 + Eta[1,0]**2
    
#def achromat_condition_print(lattice):
# #   M=np.matrix([ [1.,0.],[0.,1.] ] )
#    Eta = np.matrix([ [0.],[0.] ] )
#    for i in range(len(lattice)):
#        if lattice[i]['type']=='bend':
#            rho = lattice[i]['rho'] 
#            L = lattice[i]['length']
#            Eta=dispersion_next(Eta, lattice[i]['Mx'], L, rho)
##            Eta[0,0] = Eta[0,0] - M[1,1]*rho*(1.0-np.cos(L/rho)) - M[0,1]*np.sin(L/rho)
##            Eta[1,0] = Eta[1,0] + M[1,0]*rho*(1.0-np.cos(L/rho)) + M[0,0]*np.sin(L/rho)
#        else:
#            Eta=dispersion_next(Eta, lattice[i]['Mx'])        
#    return Eta
#    
#def achromat_condition_print2(lattice):
#    M=np.matrix([ [1.,0.],[0.,1.] ] )
#    Eta = np.matrix([ [0.],[0.] ] )
#    for i in range(len(lattice)):
#        if lattice[i]['type']=='bend':
#            rho = lattice[i]['rho'] 
#            L = lattice[i]['length']
#            Eta[0,0] = Eta[0,0] - M[1,1]*rho*(1.0-np.cos(L/rho)) - M[0,1]*np.sin(L/rho)
#            Eta[1,0] = Eta[1,0] + M[1,0]*rho*(1.0-np.cos(L/rho)) + M[0,0]*np.sin(L/rho)
#        M=lattice[i]['Mx']*M
#    Eta = M*Eta
#        
#    return Eta    
    
def SC_achromat_condition(lattice):
    M=np.matrix([ [1.,0.],[0.,1.] ] )
    EtaSC = np.matrix([ [0.],[0.] ] )
    Ltot = 0
    for i in range(len(lattice)):
        Ltot=Ltot+lattice[i]['length']
    s= -0.5*Ltot 
    for i in range(len(lattice)):
        if lattice[i]['type']=='bend':
            rho = lattice[i]['rho'] 
            L = lattice[i]['length']
            EtaSC[0,0] = EtaSC[0,0] - rho*( -M[0,1] + s*M[1,1]) +\
                                    - ( M[0,1] - (s+L)*M[1,1] )*rho*np.cos(L/rho) +\
                                    - ( (s+L)*M[0,1] + M[1,1]*rho*rho )*np.sin(L/rho) 
            EtaSC[1,0] = EtaSC[1,0] + rho*(M[0,0] - s*M[1,0]) +\
                                      ( M[0,0] - (s+L)*M[1,0] )*rho*np.cos(L/rho) +\
                                      ( (s+L)*M[0,0] + M[1,0]*rho*rho )*np.sin(L/rho) 
        M=lattice[i]['Mx']*M
        s=s+lattice[i]['length']
            
    EtaSC = M*EtaSC/Ltot*0.01  
    # 0.01 factor means deltaP/P grows about 1% trought arc length L_tot
    # if SC force is large this factor need to be increased
    return EtaSC[0,0]**2 + EtaSC[1,0]**2


#%%============================================================================
#                               design function                             
#==============================================================================
#def deltaX_at_merger(R,r,L):
#    def f(dx):
#        return R*np.cos(0.5*np.pi-L/R)-(r+dx)*np.cos(0.5*np.pi-L/r)
#    return optimize.brentq(f, 0, 0.05)  


#%%============================================================================
#                         pyNaff - rectangular window only, t=[0,T-1]
#==============================================================================
def pyNaff_full(nmode,signal,window_id=1):
    """
    tunes,amps,substracted_signals = pyNaff_full(nmode,signal)
    
    rectangular window only, 
    t=[0,T-1]
    amp = signal*np.exp(-2j*pi*tune*np.arange(T))
    """
    try:
        import scipy.optimize as opt
        from copy import deepcopy as copy
    except:
        print('util.pyNaff is not available')
        
    pi = np.pi
    T = len(signal)
    window = (1.0+np.cos(np.pi*(-1.0+2.0/(T+1.0)*np.arange(1,T+1))))**window_id
    window = window/np.sum(window)
    
    
    def getPeakInfo(signal):
        T = len(signal)
        def loss(tune):
            return -np.abs(np.sum(signal*window*np.exp(-2j*pi*tune*np.arange(T))))
        tune = np.argmax(np.abs(np.fft.fft(signal)))/T
#         result = opt.minimize(loss,tune,args=(signal),method='Nelder-Mead')
        result = opt.differential_evolution(loss,((tune-2.2/T,tune+2.2/T),),popsize=9)
        return result
    
    
    tunes = []
    amps = []
    subtracted_signals = []
    
    X = copy(signal)
    for i in range(nmode):
        result = getPeakInfo(X)
        if result.message!='Optimization terminated successfully.':
            print('Optimization failed at '+str(i+1)+'-th mode')
            break
        tunes.append(copy(result.x))
        amps.append(np.sum(X*np.exp(-2j*pi*tunes[-1]*np.arange(T)))/T)
        
        X = X - amps[-1]*np.exp(2j*pi*tunes[-1]*np.arange(T))
        subtracted_signals.append(copy(X))
    if nmode==1:
        return tunes[0],amps[0],subtracted_signals[0]
    else:
        return tunes,amps,subtracted_signals


def pyNaff(nmode,signal,window_id=1):
    tunes,amps,subtracted_signals = pyNaff_full(nmode,signal,window_id)
    return tunes,amps

#%%============================================================================
#                       (Impact lattice) closed orbit search                              
#==============================================================================


def get_tune_via_tracking(beamIn,latticeIn,pData6D_init,direction='x',nturn=1024, order=3):
    from impactIO import clearLattice
    from impactIO import writeInputFile
    from impactIO import writeParticleData
    from impactIO import run
    from impactIO import getElem
    from impactIO import readTBT
    from copy import deepcopy as copy
    
    fID = 8563
    
    beam = copy(beamIn)
    beam.nCore_y=1
    beam.nCore_z=1
    beam.n_particles=1
    beam.distribution.distribution_type = 'ReadFile'
    
    lattice = clearLattice(latticeIn)
    loop = getElem('loop')
    loop.turns = nturn
    lattice.insert(0,loop)
    TBT = getElem('TBT')
    TBT.file_id = fID
    TBT.pID_begin = 1
    TBT.pID_end = 1
    lattice.insert(1,TBT)

    data = np.zeros([1,9])
    data[0,:6]= pData6D_init
    data[0,8] = 1
    data[0,6] = beam.multi_charge.q_m[0]
    writeParticleData(data, beam.kinetic_energy, beam.mass, beam.frequency)
    
    writeInputFile(beam,lattice)
    run(order=order)    
    

    iTBT,TBT = readTBT(fID, beam.kinetic_energy, beam.mass, beam.frequency)
    if direction=='x':
        signal = TBT[:,0,0] - 1j*TBT[:,1,0]
    elif direction=='y':
        signal = TBT[:,2,0] - 1j*TBT[:,3,0]
    elif direction=='z':
        signal = TBT[:,4,0] - 1j*TBT[:,5,0]
    
    signal = signal -np.mean(signal)
    if nturn in [2**i for i in range(5,20)]:
        tune,amp,dummy = naff(1,signal,window_id=1)
    else:
        tune,amp,dummy = pyNaff(1,signal)

    return tune[0], amp[0]
    
    
        
def get_closed_orbit4D(beamIn,latticeIn,pData4D_init,delta):
    from impactIO import clearLattice
    from impactIO import writeInputFile
    from impactIO import getElem
    from copy import deepcopy as copy
    try:
        import scipy.optimize as opt
    except:
        print('util.get_closed_orbit is not available')
    
    
    fID = 5926    
    beam = copy(beamIn)
    beam.nCore_y=1
    beam.nCore_z=1
    beam.n_particles=1
    beam.distribution.distribution_type = 'ReadFile'
    gam0 = beam.kinetic_energy/beam.mass+1.0
    bet0 = np.sqrt((gam0+1.0)*(gam0-1.0))/gam0
    bg0  = gam0*bet0
    
    lattice = clearLattice(latticeIn)
    loop = getElem('loop')
    loop.turns = 1
    lattice.insert(0,loop)
    elemWrite = getElem('write_raw_ptcl')
    elemWrite.file_id = fID
    elemWrite.format_id = 2
    lattice.append(elemWrite)
    
    writeInputFile(beam,lattice)
    
    
    def _cost_closed_orbit(pData4D_init,args):
        from impactIO import writeParticleData
        from impactIO import readParticleData
        from impactIO import run
       
        fID = args[0]
        beam = args[1]
        delta= args[2]
            
        if len(args)>3:
            order=args[3]
        else:
            order=3

        data = np.zeros([1,9])
        data[0,:4]= pData4D_init
        data[0,5] = delta*(gam0+1)/gam0*beam.kinetic_energy
        data[0,8] = 1
        data[0,6] = beam.multi_charge.q_m[0]


        writeParticleData(data, beam.kinetic_energy, beam.mass, beam.frequency)
        run(order=order)
        dataOut = readParticleData(fID, beam.kinetic_energy, beam.mass, beam.frequency, format_id=2)

        diff = pData4D_init - dataOut[0,:4]
        return np.sum(diff*diff)
    
    result = opt.minimize(_cost_closed_orbit,pData4D_init,args=([fID,beam,delta]),method='Nelder-Mead',tol=1.0e-7)
    if result.message=='Optimization terminated successfully.':
        print('closed orbit optimization terminated successfully.')
        print('Closed orbit solution:',result.x)
        print('Least squre error sum:', result.fun)
        
        return result.x

    
    
import pandas as __pd
    
    
    
def getTransferMap4D(beamIn,lattice,delta=0,epsilon=[1e-8,1e-6,1e-8,1e-6],order=3):
    """
    M = getTransferMap(lattice,q,mass,ke,freq,
                       epsilon=[1e-8,1e-6,1e-8,1e-6,1e-7,1.0],
                       fname='test.in' )
    get linear transfer map (without space-charge)  by tracking 6 particles
    whose initial phase-space perturbation given by epsilon
    input
        beamIn  = impact beam class
        lattice = (dict) lattice dictionary whose transvermap to be determined
        epsilon = 6 dimension array of perturbation for 
                  x,px,y,py, z*360/v/freq, E  in unit of 
                  [m],[rad],[m],[rad],[deg],[MeV]
                  default : epsilon = [1e-e-8,1e-6,1e-8,1e-6,1e-7,1.0]
    """
    from impactIO import writeInputFile
    from impactIO import writeParticleData
    from impactIO import readParticleData
    from impactIO import run
    from impactIO import getElem
    from copy import deepcopy as copy
    
    
    
    if delta!=0:
        closed_orbit = get_closed_orbit4D(beamIn,lattice,np.zeros(4),delta)
    else:
        closed_orbit = np.zeros(4)
    
    beam = copy(beamIn)
    beam.nCore_y=1
    beam.nCore_z=1
    beam.n_particles=6
    beam.distribution.distribution_type = 'ReadFile'
    gam0 = beam.kinetic_energy/beam.mass+1.0
    
    line = lattice.copy()
    line = [item for item in line if not item.type == 'write_raw_ptcl']    
    line = [item for item in line if not item.type == 'loop']   
    loop = getElem('loop')
    loop.turns = 1
    line.insert(0,loop)
    elemWrite = getElem('write_raw_ptcl')
    elemWrite.file_id = 5926
    elemWrite.format_id = 2
    line.append(elemWrite)
    writeInputFile(beam,line)
    
    data = np.zeros([4,9])
    for i in range(4):
        data[i,:4] = closed_orbit[:]
    for i in range(4):
        data[i,i] = data[i,i] + epsilon[i]
        data[i,8] = i+1
    data[:,5] = delta*(gam0+1)/gam0*beam.kinetic_energy
    data[:,6] = beam.multi_charge.q_m[0]
    writeParticleData(data,beam.kinetic_energy, beam.mass, beam.frequency)
    run(order=order)
    
    dataOut = readParticleData(5926, beam.kinetic_energy, beam.mass, beam.frequency, format_id=2)[:,:4]
    #os.system('rm fort.'+str(fileID))
    m,n = dataOut.shape
    M = np.zeros([4,4])
    if m<4:
        print('particle lost observed. too large inital perturbation')
    else:
        for i in range(4):
            M[:,i] = (dataOut[i,:]-closed_orbit)/epsilon[i]
    return __pd.DataFrame(M)


    
def getTransferMap6D(beamIn,lattice,epsilon=[1e-8,1e-6,1e-8,1e-6,1e-7,1.0],order=3):
    """
    M = getTransferMap(lattice,q,mass,ke,freq,
                       epsilon=[1e-8,1e-6,1e-8,1e-6,1e-7,1.0],
                       fname='test.in' )
    get linear transfer map (without space-charge)  by tracking 6 particles
    whose initial phase-space perturbation given by epsilon
    input
        beamIn  = impact beam class
        lattice = (dict) lattice dictionary whose transvermap to be determined
        epsilon = 6 dimension array of perturbation for 
                  x,px,y,py, z*360/v/freq, E  in unit of 
                  [m],[rad],[m],[rad],[deg],[MeV]
                  default : epsilon = [1e-e-8,1e-6,1e-8,1e-6,1e-7,1.0]
    """
    from impactIO import writeInputFile
    from impactIO import writeParticleData
    from impactIO import readParticleData
    from impactIO import run
    from impactIO import getElem
    from copy import deepcopy as copy
    
    
    beam = copy(beamIn)
    beam.nCore_y=1
    beam.nCore_z=1
    beam.n_particles=6
    beam.distribution.distribution_type = 'ReadFile'
    gam0 = beam.kinetic_energy/beam.mass+1.0
    bet0 = np.sqrt((gam0+1.0)*(gam0-1.0))/gam0
    bg0  = gam0*bet0
    
    line = lattice.copy()
    line = [item for item in line if not item.type == 'write_raw_ptcl']    
    line = [item for item in line if not item.type == 'loop']   
    loop = getElem('loop')
    loop.turns = 1
    line.insert(0,loop)
    elemWrite = getElem('write_raw_ptcl')
    elemWrite.file_id = 5926
    elemWrite.format_id = 2
    line.append(elemWrite)
    writeInputFile(beam,line)
    
    data = np.zeros([6,9])
    for i in range(6):
        data[i,i] = epsilon[i]
        data[i,8] = i+1
    data[:,6] = beam.multi_charge.q_m[0]
    writeParticleData(data,beam.kinetic_energy, beam.mass, beam.frequency)
    run(order=order)
    
    dataOut = readParticleData(5926, beam.kinetic_energy, beam.mass, beam.frequency, format_id=2)[:,:6]
    #os.system('rm fort.'+str(fileID))
    m,n = dataOut.shape
    M = np.zeros([6,6])
    if m<6:
        print('particle lost observed. too large inital perturbation')
    else:
        for i in range(6):
            M[:,i] = dataOut[i,:]/epsilon[i]
    return __pd.DataFrame(M)