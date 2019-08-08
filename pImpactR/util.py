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
