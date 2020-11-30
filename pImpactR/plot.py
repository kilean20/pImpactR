# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 12:48:25 2016

@author: kilean
"""

import numpy as np
try:
    import matplotlib.pyplot as plt
except:
    print('matplotlib not found. plot module is disabled')
import os
#%%############################################################################
##            internel class and function to layout lattice plot             ##
###############################################################################
class _element:
    def __init__(self, lattice, elem_id):
        self.number=0
        self.elem_id=elem_id
        
        for i in range(len(lattice)):
            if lattice[i,1]==self.elem_id:
                self.number=self.number+1;
        self.positions=np.zeros([4*self.number+2,2])
        count=0
        z=0.0
        for i in range(0,len(lattice)):
            if lattice[i,1]==self.elem_id:
                self.positions[4*count+1,0]=z
                self.positions[4*count+2,0]=z+1E-15;
                self.positions[4*count+3,0]=z+lattice[i,0]-1E-15;
                self.positions[4*count+4,0]=z+lattice[i,0];
                if lattice[i,2]>0:
                    self.positions[4*count+2,1]=1;
                    self.positions[4*count+3,1]=1;
                else :
                    self.positions[4*count+2,1]=-1;
                    self.positions[4*count+3,1]=-1;                    
                count=count+1
            z=z+lattice[i,0]
        self.positions[0,0]=0.0
        self.positions[-1,0]=z               

             
def _lattice(h, lattice_info, fileDir=''):
    
    lattice_offset = lattice_info['offset']
    lattice_scale = lattice_info['scale']
    if 'file' in lattice_info:
        lattice_file = lattice_info['file']
    else:
        lattice_file = 'test.in'
    
    lattice=np.loadtxt(fileDir+lattice_file,
                       comments='!',
                       skiprows=12,
                       usecols=(0,3,4))
                                                  
    RF=_element(lattice,104)
    Bend=_element(lattice,4)
    Quad=_element(lattice,1)
    
    
    h.fill_between(RF.positions[:,0], lattice_offset, 
                   lattice_offset+lattice_scale*RF.positions[:,1], alpha=0.2, 
                   facecolor='black', color="none", hatch="X")
    h.fill_between(RF.positions[:,0], lattice_offset, 
                   lattice_offset-lattice_scale*RF.positions[:,1], alpha=0.2, 
                   facecolor='black', color="none", hatch="X")
    h.fill_between(Quad.positions[:,0], lattice_offset, 
                   lattice_offset+lattice_scale*Quad.positions[:,1], alpha=0.6, 
                   facecolor='black', linewidth=0.2)
    h.fill_between(Bend.positions[:,0], lattice_offset, 
                   lattice_offset+lattice_scale*Bend.positions[:,1], alpha=0.1, 
                   facecolor='black', color="none", hatch="o")
    L=0.0
    for elem in lattice:
      L=L+elem[0]
    h.plot([0,L],[lattice_offset,lattice_offset],'k-',linewidth=0.5)
def _getLatticeInfo(lattice_info=None, minX=0.0, maxX=1.0):
    if lattice_info==None:
        return {'file':'test.in','offset':None, 'scale':None}
  
    lattice_info['offset'] = 0.5*(minX+maxX)
    lattice_info['scale'] = 0.2*( 0.5*(maxX-minX) )
    
    
#%%############################################################################
###############################################################################
###                                  RMS plots                              ###
###############################################################################
############################################################################### 

#%%#####  internel functions  ######
# xy envelope
def _x_envelope(h, fileDir='',flag_maxRadius=False,
                plotRange=None, lattice_info=None):
    '''                
    h = matplotlib handle
    '''
    X=np.loadtxt(fileDir+'fort.24')
    X_centroid=X[:,1]*1E3
    X_rms=X[:,2]*1E3
    
    Y=np.loadtxt(fileDir+'fort.25')
    Y_centroid=Y[:,1]*1E3
    Y_rms=Y[:,2]*1E3
    
    S=X[:,0]
    s_max=max(S)
    s_min=min(S) 
    
    h.plot(S,X_centroid+X_rms,'b',S,X_centroid-X_rms,'b')
    h.plot(S,Y_centroid+Y_rms,'r',S,Y_centroid-Y_rms,'r')
    h.set_ylabel(r'$\mathsf{\sigma_{x,y} \, (mm)}$', fontsize=12)    
    
    if lattice_info != None:
        if lattice_info['scale']==None:
            minX = min( min(X_centroid -X_rms),  min(Y_centroid -Y_rms) )
            maxX = max( max(X_centroid +X_rms),  max(Y_centroid +Y_rms) )
            _getLatticeInfo(lattice_info, minX, maxX)
        _lattice(h, lattice_info, fileDir=fileDir) 
    
    if flag_maxRadius:
        Rmax=1E3*np.loadtxt(fileDir+'fort.18')[:,-1]
        h.plot(S,-Rmax,'g',S,Rmax,'g')
        
    if plotRange==None:
        h.set_xlim([s_min,s_max])
    else:
        if 's' in plotRange:
            h.set_xlim(plotRange['s'])        
        if 'x' in plotRange:        
            h.set_ylim(plotRange['x'])

# p_xy envelope
def _px_envelope(h, fileDir='', plotRange=None, lattice_info=None):
    '''
    in/out : h = matplotlib handle
    '''
    X=np.loadtxt(fileDir+'fort.24')
    Y=np.loadtxt(fileDir+'fort.25')
    
    Px_centroid=X[:,3]*1E3
    Px_rms=X[:,4]*1E3
    Py_centroid=Y[:,3]*1E3
    Py_rms=Y[:,4]*1E3
    
    S=X[:,0]
    s_max=max(S)
    s_min=min(S)   
    
    h.plot(S,Px_centroid+Px_rms,'b',S,Px_centroid-Px_rms,'b')
    h.plot(S,Py_centroid+Py_rms,'r',S,Py_centroid-Py_rms,'r')
    h.set_ylabel(r'$\mathsf{\sigma_{px,py}\, (mrad)}$', fontsize=12)    
    
    if lattice_info != None:
        if lattice_info['scale']==None:
            minX = min( min(Px_centroid -Px_rms),  min(Py_centroid -Py_rms) )
            maxX = max( max(Px_centroid +Px_rms),  max(Py_centroid +Py_rms) )
            _getLatticeInfo(lattice_info, minX, maxX)
        _lattice(h, lattice_info, fileDir=fileDir) 
    
    if plotRange==None:
        h.set_xlim([s_min,s_max])
    else:
        if 's' in plotRange:
            h.set_xlim(plotRange['s'])        
        if 'px' in plotRange:        
            h.set_ylim(plotRange['px'])

def _x_emittance(h, halo=None, fileDir='', plotRange=None, lattice_info=None):
    '''                
    h = matplotlib handle
    '''
    X=np.loadtxt(fileDir+'fort.24')
    X_emittance=X[:,6]*1E6
    
    Y=np.loadtxt(fileDir+'fort.25')
    Y_emittance=Y[:,6]*1E6
    
    S=X[:,0]
    s_max=max(S)
    s_min=min(S)   
    
    h.plot(S,X_emittance,'b')
    h.plot(S,Y_emittance,'r')
    h.set_ylabel(r'$\mathsf{norm.\, \epsilon \, (mm \, mrad)}$', fontsize=11)   
    
    if lattice_info != None:
        if lattice_info['scale']==None:
            minX = min( min(X_emittance),  min(Y_emittance) )
            maxX = max( max(X_emittance),  max(Y_emittance) )
            _getLatticeInfo(lattice_info, minX, maxX)
        _lattice(h, lattice_info, fileDir=fileDir)    
    
    if isinstance(halo,int) and X.shape[1] == 10:
        h2 = h.twinx()
        if halo==99:
            X_halo=X[:,9]*1E6
            Y_halo=Y[:,9]*1E6
        elif halo==90:
            X_halo=X[:,7]*1E6
            Y_halo=Y[:,7]*1E6
        else :
            X_halo=X[:,8]*1E6
            Y_halo=Y[:,8]*1E6
        
        h2.plot(S,1-X_halo/X_emittance,'b--',alpha=0.5)
        h2.plot(S,1-Y_halo/Y_emittance,'r--',alpha=0.5)
        h2.set_ylabel(r'$\mathsf{ \Delta \epsilon / \epsilon}$', fontsize=12)           
    elif isinstance(halo,int):
        print('90,95,99% rms emittance is not calculate. IMPACT standard output flagg must be set to 2')

    if plotRange==None:
        h.set_xlim([s_min,s_max])
    else:
        if 's' in plotRange:
            h.set_xlim(plotRange['s'])
            if isinstance(halo,int) and X.shape[1] == 10:
                h2.set_xlim(plotRange['s'])  
        if 'ex' in plotRange:        
            h.set_ylim(plotRange['ex']) 
        if isinstance(halo,int) and X.shape[1] == 10 and 'Delta_ex' in plotRange:
            h2.set_ylim(plotRange['Delta_ex'])              
    

#%% Plot longitudinal RMS statistics
def _z_envelope(h, fileDir='', plotRange=None, lattice_info=None):
    '''                
    h = matplotlib handle
    '''
    Z=np.loadtxt(fileDir+'fort.26')
    Z_centroid=Z[:,1]
    Z_rms=Z[:,2]
    
    S=Z[:,0]
    s_max=max(S)
    s_min=min(S) 
    
    h.plot(S,Z_centroid+Z_rms,'b',S,Z_centroid-Z_rms,'b')
    h.set_ylabel(r'$\mathsf{\sigma_z \, (deg)}$', fontsize=12)
    
    if lattice_info != None:
        if lattice_info['scale']==None:
            minX = min(Z_centroid -Z_rms)
            maxX = max(Z_centroid +Z_rms)
            _getLatticeInfo(lattice_info, minX, maxX)
        _lattice(h, lattice_info, fileDir=fileDir) 

    if plotRange==None:
        h.set_xlim([s_min,s_max])
    else:
        if 's' in plotRange:
            h.set_xlim(plotRange['s'])        
        if 'z' in plotRange:        
            h.set_ylim(plotRange['z'])
            
    
def _pz_envelope(h, fileDir='', plotRange=None, lattice_info=None):
    '''                
    h = matplotlib handle
    '''
    Z=np.loadtxt(fileDir+'fort.26')
    Pz_centroid=Z[:,3]
    Pz_rms=Z[:,4]
    S=Z[:,0]
    s_max=max(S)
    s_min=min(S)   
    
    h.plot(S,Pz_centroid+Pz_rms,'b',S,Pz_centroid-Pz_rms,'b')
    h.set_ylabel(r'$\mathsf{\sigma_E \, (MeV)}$', fontsize=12)
    
    if lattice_info != None:
        if lattice_info['scale']==None:
            minX = min(Pz_centroid -Pz_rms)
            maxX = max(Pz_centroid +Pz_rms)
            _getLatticeInfo(lattice_info, minX, maxX)
        _lattice(h, lattice_info, fileDir=fileDir) 
    
    if plotRange==None:
        h.set_xlim([s_min,s_max])
    else:
        if 's' in plotRange:
            h.set_xlim(plotRange['s'])        
        if 'px' in plotRange:        
            h.set_ylim(plotRange['px'])         

     
def _z_emittance(h, halo=None, fileDir='', plotRange=None, lattice_info=None):
    '''                
    h = matplotlib handle
    '''
    Z=np.loadtxt(fileDir+'fort.26')
    Z_emittance=Z[:,6]   
    S=Z[:,0]
    s_max=max(S)
    s_min=min(S)   
    
    h.plot(S,Z_emittance,'b')
    h.set_ylabel(r'$\mathsf{\epsilon \, (deg\,MeV)}$', fontsize=11)
    
    if lattice_info != None:
        if lattice_info['scale']==None:
            minX = min(Z_emittance)
            maxX = max(Z_emittance)
            _getLatticeInfo(lattice_info, minX, maxX)
        _lattice(h, lattice_info, fileDir=fileDir) 
     
    
    if isinstance(halo,int) and Z.shape[1] == 10:
        h2 = h.twinx()
        if halo==99:
            Z_halo=Z[:,9]
        elif halo==90:
            Z_halo=Z[:,7]
        else :
            Z_halo=Z[:,8]           
        h2.plot(S,1-Z_halo/Z_emittance,'b--',alpha=0.5)
        h2.set_ylabel(r'$\mathsf{ \Delta \epsilon / \epsilon}$', fontsize=12)  
    elif isinstance(halo,int):
        print('90,95,99% rms emittance is not calculate. IMPACT standard output flagg must be set to 2')
        
    if plotRange==None:
        h.set_xlim([s_min,s_max])
    else:
        if 's' in plotRange:
            h.set_xlim(plotRange['s'])   
            if isinstance(halo,int) and Z.shape[1] == 10:
                h2.set_xlim(plotRange['s'])  
        if 'ez' in plotRange:        
            h.set_ylim(plotRange['ez'])  
        if isinstance(halo,int) and Z.shape[1] == 10 and 'Delta_ez' in plotRange:
            h2.set_ylim(plotRange['Delta_ez'])   

    
    
#%% Plot RMS statistics
def _getLatticeInfo_4rmsPlot(x=True, px=True, ex=False, z=True, pz=True, ez=False):
    """
    lattice_info = _getLatticeInfo_4rmsPlot(x=True, px=True, ex=False, z=True, pz=True, ez=False)
    """
    lattice_info = {}
    if x:
        lattice_info['x'] = {'file':'test.in','offset':None, 'scale':None}
    if px:
        lattice_info['px'] = {'file':'test.in','offset':None, 'scale':None}
    if ex:
        lattice_info['ex'] = {'file':'test.in','offset':None, 'scale':None}
    if z:
        lattice_info['z'] = {'file':'test.in','offset':None, 'scale':None}
    if pz:
        lattice_info['pz'] = {'file':'test.in','offset':None, 'scale':None}
    if ez:
        lattice_info['ez'] = {'file':'test.in','offset':None, 'scale':None}
    return lattice_info    
    
def rms(savefileID=0, fileDir='', lattice_info=_getLatticeInfo_4rmsPlot(), flag_maxRadius=False,
        plot_momentum=True, halo=None, plotRange=None):
    """
    figX,figZ = rms(savefileID=0, fileDir='', 
                    lattice_info=_getLatticeInfo_4rmsPlot(), 
                    flag_maxRadius=False, plot_momentum=True,
                    halo=None, plotRange=None)
    plot rms beam profile evolution. IMPACT must be ran priori.
    input 
        savefileID : (int) fileID for plot save
        plot_momentum : (bool) plot px,py,pz or not
        flag_maxRadius : (bool) plot max amplitude out of all particles
        plotRange : (matplotlib plot range)
        halo : (int) any of 90,95,99 indicating 90,95,99% rms quantities
               (priori : non-standard IMPACT ouput flag)
        lattice_info : for advanced user
                       infer the source of _getLatticeInfo_4rmsPlot() 
    output
      figX : matplotlib figure for x,y RMS info
      figY : matplotlib figure for z RMS info
    output to file
      fileDir+'x'+str(savefileID)+'.png'
      fileDir+'z'+str(savefileID)+'.png'
    """            
    plt.rcParams.update({'font.size': 9})  
    # transverse
    figX = plt.figure() 
    
    if plot_momentum:
        plotid = 311
    else:
        plotid = 211
    sub1=figX.add_subplot(plotid)
    if 'x' in lattice_info:
        info = lattice_info['x']
    else :
        info = None
    _x_envelope(sub1,fileDir=fileDir, flag_maxRadius=flag_maxRadius,
                 plotRange=plotRange, lattice_info=info)


    if plot_momentum:
        plotid=plotid+1
        if 'px' in lattice_info:
            info = lattice_info['px']
        else :
            info = None
        sub2=figX.add_subplot(plotid)
        _px_envelope(sub2,fileDir=fileDir, plotRange=plotRange, lattice_info=info)


    plotid=plotid+1
    sub3=figX.add_subplot(plotid)
    if 'ex' in lattice_info:
        info = lattice_info['ex']
    else :
        info = None
    _x_emittance(sub3,halo=halo,fileDir=fileDir, plotRange=plotRange, lattice_info=info)
    sub3.set_xlabel('z (m)')
    
    
    plt.tight_layout() 
    plt.savefig(fileDir+'x'+str(savefileID)+'.png', dpi=240)
    plt.close()
    
    # longitudinal
    figZ = plt.figure() 
    if plot_momentum:
        plotid = 311
    else:
        plotid = 211
    sub1=figZ.add_subplot(plotid)
    if 'z' in lattice_info:
        info = lattice_info['z']
    else :
        info = None
    _z_envelope(sub1,fileDir=fileDir, 
                 plotRange=plotRange, lattice_info=info)


    if plot_momentum:
        plotid=plotid+1
        sub2=figZ.add_subplot(plotid)
        if 'pz' in lattice_info:
            info = lattice_info['pz']
        else :
            info = None
        _pz_envelope(sub2,fileDir=fileDir, plotRange=plotRange, lattice_info=info)


    plotid=plotid+1
    sub3=figZ.add_subplot(plotid)
    if 'ez' in lattice_info:
        info = lattice_info['ez']
    else :
        info = None    
    _z_emittance(sub3,halo=halo,fileDir=fileDir, plotRange=plotRange, lattice_info=info)
    sub3.set_xlabel('z (m)')


    plt.tight_layout() 
    plt.savefig(fileDir+'z'+str(savefileID)+'.png', dpi=240)
    plt.close()
    return figX,figZ
  
#%%############################################################################
###############################################################################
###                            max amplitude plots                          ###
###############################################################################
############################################################################### 
def _getLatticeInfo_4maxPlot(x=True, px=True, z=True, pz=True):
    lattice_info = {}
    if x:
        lattice_info['x'] = {'file':'test.in','offset':None, 'scale':None}
    if px:
        lattice_info['px'] = {'file':'test.in','offset':None, 'scale':None}
    if z:
        lattice_info['z'] = {'file':'test.in','offset':None, 'scale':None}
    if pz:
        lattice_info['pz'] = {'file':'test.in','offset':None, 'scale':None}
    return lattice_info    

def maxAmplitude(savefileID=0, fileDir='', 
        plot_momentum=True, flag_radius=True,
        plotRange=None, lattice_info=_getLatticeInfo_4maxPlot() ):
    
    plt.rcParams.update({'font.size': 9})            

    Data=np.loadtxt(fileDir+'fort.27')*1E3
    S  = Data[:,0]*1E-3
    X  = Data[:,1]
    Px = Data[:,2]
    Y  = Data[:,3]
    Py = Data[:,4]
    R  = np.loadtxt(fileDir+'fort.18')[:,-1]*1E3
    Z  = Data[:,5]*1E-3
    Pz = Data[:,6] 
    s_max=max(S)
    s_min=min(S)   
    
    plotid=111
    fig = plt.figure() 
    sub1=fig.add_subplot(plotid)
    if flag_radius:
        sub1.plot(S,R,'b')
        sub1.set_ylabel(r'$\mathsf{ R_{max} \, (mm)}$', fontsize=12)
    else:
        sub1.plot(S,X,'b',S,Y,'b--',alpha=0.7)
        sub1.set_ylabel(r'$\mathsf{ x,y_{max} \, (mm)}$', fontsize=12)
    sub1.set_xlabel('z (m)')
        
    if lattice_info['x'] != None:
        if lattice_info['x']['scale']==None:
            minX = min(min(X),min(Y)) 
            maxX = max(max(X),max(Y)) 
            _getLatticeInfo(lattice_info['x'], minX, maxX)
        _lattice(sub1, lattice_info['x'], fileDir=fileDir) 
        
    
    sub2 = sub1.twinx()  
    sub2.plot(S,Z,'r')
    sub2.set_ylabel(r'$\mathsf{ z_{max} \, (deg)}$', fontsize=12)       
        
    if plotRange==None:
        sub1.set_xlim([s_min,s_max])
    else:
        if 's' in plotRange:
            sub1.set_xlim(plotRange['s'])   
            sub2.set_xlim(plotRange['s']) 
        if 'x' in plotRange:   
            sub1.set_ylim(plotRange['x'])   
        if 'z' in plotRange:     
            sub2.set_ylim(plotRange['z'])
        
    plt.tight_layout() 
    plt.savefig(fileDir+'max_radius'+str(savefileID)+'.png', dpi=240)



    fig = plt.figure() 
    sub1=fig.add_subplot(plotid)
    sub1.plot(S,Pz,'b')
    sub1.set_ylabel(r'$\mathsf{ \Delta E_{max} \, (MeV)}$', fontsize=12)
    sub1.set_xlabel('z (m)')
    
    if lattice_info['pz'] != None:
        if lattice_info['pz']['scale']==None:
            minX = min(Pz)
            maxX = max(Pz)  
            _getLatticeInfo(lattice_info['pz'], minX, maxX)
            _lattice(sub1, lattice_info['pz'], fileDir=fileDir) 
            
    if plotRange==None:
        sub1.set_xlim([s_min,s_max])
    else:
        if 's' in plotRange:
            sub2.set_xlim(plotRange['s'])        
        if 'pz' in plotRange:        
            sub2.set_ylim(plotRange['pz'])
 
    if plot_momentum:
        sub2 = sub1.twinx()  
        sub2.plot(S,Px,'r',S,Py,'r--',alpha=0.7)
        sub2.set_ylabel(r'$\mathsf{ p_{x,y,max} \, (mrad)}$', fontsize=12)
        if plotRange !=None and 'px' in plotRange:
            sub2.set_ylim(plotRange['px'])        

    plt.tight_layout() 
    plt.savefig(fileDir+'max_energy'+str(savefileID)+'.png', dpi=240)
    plt.close()        

        



    
    
#%%############################################################################
###############################################################################
###                               particle loss                             ###
###############################################################################
############################################################################### 
from impactIO import readInputFile
def powerLoss(savefileID=0, fileDir='', lattice_info=None, plotRange=None):
            
    
    dummy = np.loadtxt(fileDir+'fort.18')
    S = dummy [:,0]
    Energy = dummy[:,3]    
    Particles = np.loadtxt(fileDir+'fort.32')[:,1]
    n_total = Particles[0]
     
    beam,lattice=readInputFile('test.in')
    PowerLoss = (1.0-Particles/n_total)*Energy*beam['current'] #assume all buckets are filled (Watt)
    
    dP = np.gradient(PowerLoss)
    dS = np.gradient(S)
    PowerLossPerMeter = dP/dS
    
    plt.rcParams.update({'font.size': 9})  
    fig = plt.figure() 
    sub1=fig.add_subplot(111)
    sub1.plot(S,PowerLossPerMeter,'r')
    ###
    if lattice_info != None:
        if lattice_info['scale']==None:
            minX = min(PowerLossPerMeter)
            maxX = max(PowerLossPerMeter) 
            _getLatticeInfo(lattice_info, minX, maxX)
        _lattice(sub1, lattice_info, fileDir=fileDir) 

    sub1.set_xlabel('z (m)')
    sub1.set_ylabel('power loss (W/m)')

    plt.tight_layout() 
    plt.savefig(fileDir+'loss'+str(savefileID)+'.png', dpi=240)
    plt.close()


#%%############################################################################
###############################################################################
###                              phase-space plots                          ###
###############################################################################
############################################################################### 
#%% plot poincare of sliced bunch
from impactIO import readParticleDataSliced
def phase_space(fileID, ke, mass, freq, zSliced=True, nSlice=1,
                plotRange = None, fileDir='', saveDir='', showAll=False):
    """
    Fig = phase_space(fileID, ke, mass, freq, zSliced=True, nSlice=1,
                      plotRange = None, fileDir='', saveDir='', 
                      showAll=False):
                
    phase-space plot of sliced particle data. 
    Density plot imposed if Scipy is installed.
    input : 
        fileID : (int) fileID of IMPACT particle output
        ke  : (real) reference kinetic energy 
        mass : (real) particle mass
        freq : (real) reference frequency used to define 
                      longitudinal coordinate
        zSliced : (bool) longitudinal / energy slice
        nSlice : (int) number of slices
        ...
    output :
       Fig : (list) matplotlib figures of each slice
    """
    pData = readParticleDataSliced(nSlice, fileID, ke, mass, freq, zSliced)
    Fig = []
    for i in range(nSlice):
        data = pData[i]
        n_particles = len(data)
        if n_particles < 10:
            continue
        
        sigmax=np.std(data[:,0])
        sigmapx=np.std(data[:,1])
        dx=data[:,0]-np.mean(data[:,0])
        dpx=data[:,1]-np.mean(data[:,1])
        mux=-np.mean(dx*dpx)
        data[:,6]= np.sqrt( dx*dx*(sigmapx*sigmapx) 
                            +2.0*mux*dx*dpx+
                            +dpx*dpx*(sigmax*sigmax) )
                
        dy=data[:,2]-np.mean(data[:,2])
        dpy=data[:,3]-np.mean(data[:,3])
        sigmay=np.std(data[:,2])
        sigmapy=np.std(data[:,3])
        muy=-np.mean(dy*dpy)                 
        data[:,7]= np.sqrt( dy*dy*(sigmapy*sigmapy) 
                            +2.0*muy*dy*dpy+
                            +dpy*dpy*(sigmay*sigmay) )
                          
        muxy=-np.mean(dx*dy);         
        data[:,8]= np.sqrt( dx*dx*(sigmay*sigmay) 
                            +2.0*muxy*dx*dy+
                            +dy*dy*(sigmax*sigmax) )
        
        if showAll :
            row_col_i = 141
            Fig.append(plt.figure(figsize=(9.6, 2.4)))
        else:
            row_col_i = 131   
            Fig.append(plt.figure(figsize=(7.2, 2.4)))
        
        plt.rcParams.update({'font.size': 9})
        ax=Fig[-1].add_subplot(row_col_i)
        ax.scatter(1E3*data[:,0], 1E3*data[:,1],
                   c=-data[:,6],marker='.',edgecolors='none' )
        ax.grid()          
        ax.set_xlabel(r"$\mathsf{x}$"+' (mm)',fontsize=10) 
        ax.set_ylabel(r"$\mathsf{p_x}$"+' (mrad)',fontsize=10) 
        
        
        ay=Fig[-1].add_subplot(row_col_i+1)
        ay.scatter(1E3*data[:,2], 1E3*data[:,3],
                   c=-data[:,7],marker='.',edgecolors='none' )
        ay.grid() 
        ay.set_xlabel(r"$\mathsf{y}$"+' (mm)',fontsize=10) 
        ay.set_ylabel(r"$\mathsf{p_y}$"+' (mrad)',fontsize=10) 

        
        axy=Fig[-1].add_subplot(row_col_i+2)
        axy.scatter(1E3*data[:,0], 1E3*data[:,2],
                   c=-data[:,8], marker='.', edgecolors='none' )
        axy.grid()
        axy.set_xlabel(r"$\mathsf{x}$"+' (mm)',fontsize=10) 
        axy.set_ylabel(r"$\mathsf{y}$"+' (mm)',fontsize=10) 

        if showAll :
            sigmaz=np.std(data[:,4])
            sigmapz=np.std(data[:,5])
            dz=data[:,4]-np.mean(data[:,4])
            dpz=data[:,5]-np.mean(data[:,5])
            muz=-np.mean(dz*dpz); 

            data[:,8]= np.sqrt( dz*dz*(sigmapz*sigmapz) 
                                +2.0*muz*dz*dpz+
                                +dpz*dpz*(sigmaz*sigmaz) )                  
                      
            az=Fig[-1].add_subplot(row_col_i+3)
            az.scatter(data[:,4], data[:,5],
                       c=-data[:,8], marker='.', edgecolors='none' )
            az.grid()
            az.set_xlabel(r"$\mathsf{\phi}$"+' (degree)',fontsize=10) 
            az.set_ylabel(r"$\mathsf{\Delta E}$"+' (MeV)',fontsize=10) 

        if plotRange != None:
            ax.set_xlim(plotRange[0])
            ax.set_ylim(plotRange[1])
            ay.set_xlim(plotRange[2])
            ay.set_ylim(plotRange[3])
            axy.set_xlim(plotRange[4])
            axy.set_ylim(plotRange[5])
            
        plt.tight_layout()
        
        if zSliced==True:
            if nSlice !=1:
                Fig[-1].suptitle( str(n_particles)+" particles of bunch slice at "+
                              "{0:.3f}".format(np.mean(data[:,4]))+" degree", 
                              y=1.05, fontsize=13)
            if saveDir == '':
                plt.savefig('phase_space'+str(fileID)+'_z'+str(i)+'.png', dpi=240)
            else:        
                if not os.path.exists(saveDir):
                    os.makedirs(saveDir)
                plt.savefig('./'+saveDir+'/phase_space'+str(fileID)+'_z'+str(i)+'.png', dpi=240)                          
        else:
            if nSlice !=1:
                Fig[-1].suptitle( str(n_particles)+" particles of bunch slice at "+
                              "{0:.3f}".format(np.mean(data[:,5]))+"MeV", 
                               y=1.05, fontsize=13)
            if saveDir == '':
                plt.savefig('phase_space'+str(fileID)+'_E'+str(i)+'.png', dpi=240)
            else:        
                if not os.path.exists(saveDir):
                    os.makedirs(saveDir)
                plt.savefig('./'+saveDir+'/phase_space'+str(fileID)+'_E'+str(i)+'.png', dpi=240)
    return Fig


#%% density plot
try:
    from scipy import stats
    def density(Xin,Yin, samplePeriod=1, xlim=None, ylim=None, xlabel=None, ylabel=None, mksize=4, pltHandle=plt):
        index_regular = ~(np.isnan(Xin)+np.isnan(Yin))
        X = Xin[index_regular]
        Y = Yin[index_regular]
        Xk=X[0::samplePeriod];Yk=Y[0::samplePeriod]
        kernel = stats.gaussian_kde([Xk,Yk])
        cData = kernel.evaluate([X,Y])
        scatter = pltHandle.scatter(X,Y, c=cData, s=mksize, lw = 0)
        pltHandle.colorbar(scatter)
#         if xlim==None:
#             xlim = [min(X),max(X)]
#         pltHandle.xlim(xlim[0],xlim[1])
#         if ylim==None:   
#             ylim = [min(Y),max(Y)]
#         pltHandle.ylim(ylim[0],ylim[1])
#         if xlabel!=None:
#             pltHandle.xlabel(xlabel)
#         if ylabel!=None:
#             pltHandle.ylabel(ylabel)
        return scatter

    def contour(x,y,levels=0,edge=True, pltHandle=plt):
      """Draw a two dimensional kernel density plot.
      You can specify either a figure or an axis to draw on.

      Parameters:
      -----------
      x: 1D numpy array
      y: 1D numpy array of same length with x

      Returns:
      --------
      fig, ax: matplotlib figure and axis objects
      """

      x_min = x.min()
      x_max = x.max()
      y_min = y.min()
      y_max = y.max()
      X, Y = np.mgrid[x_min:x_max:200j, y_min:y_max:200j]
      positions = np.vstack([X.ravel(), Y.ravel()])
      values = np.vstack([x, y])
      import scipy.stats as stats
      kernel = stats.gaussian_kde(values)
      Z = np.reshape(kernel(positions).T, X.shape)
      if edge:
        maxLevel = np.max(kernel)
        levels = np.linspace(0.01*maxLevel,maxLevel,10)
        return pltHandle.contour(Z, extent=[x_min, x_max, y_min, y_max],levels=levels)
      if np.all(levels==0):
        return pltHandle.contour(Z, extent=[x_min, x_max, y_min, y_max])
      else:
        return pltHandle.contour(Z, extent=[x_min, x_max, y_min, y_max],levels=levels)
            
except:
    print('scipy not found. plot.density module is disabled')

        
        
try:
    import matplotlib.colors as mplc
    from matplotlib.colors import colorConverter as cC
except:
    print('matplotlib not found. plot module is disabled')
        


def _to_rgb(c):
    """
    Convert color *c* to a numpy array of *RGB* handling exeption
    Parameters
    ----------
    c: Matplotlib color
        same as *color* in *colorAlpha_to_rgb*
    output
    ------
    rgbs: list of numpy array
        list of c converted to *RGB* array
    """

    if(getattr(c, '__iter__', False) == False):  #if1: if c is a single element (number of string)
        rgbs = [np.array(cC.to_rgb(c)),]  #list with 1 RGB numpy array

    else:  #if1, else: if is more that one element

        try:   #try1: check if c is numberic or not
            np.array(c) + 1

        except (TypeError, ValueError):  #try1: if not numerics is not (only) RGB or RGBA colors
            #convert the list/tuble/array of colors into a list of numpy arrays of RGB
            rgbs = [np.array( cC.to_rgb(i)) for i in c]

        except Exception as e:  #try1: if any other exception raised
            print("Unexpected error: {}".format(e))
            raise e #raise it

        else:  #try1: if the colors are all numberics

            arrc = np.array(c)  #convert c to a numpy array
            arrcsh = arrc.shape  #shape of the array 

            if len(arrcsh)==1:  #if2: if 1D array given 
                if(arrcsh[0]==3 or arrcsh[0]==4):  #if3: if RGB or RBGA
                    rgbs = [np.array(cC.to_rgb(c)),]  #list with 1 RGB numpy array
                else:   #if3, else: the color cannot be RBG or RGBA
                    raise ValueError('Invalid rgb arg "{}"'.format(c))
                #end if3
            elif len(arrcsh)==2:  #if2, else: if 2D array
                if(arrcsh[1]==3 or arrcsh[1]==4):  #if4: if RGB or RBGA
                    rgbs = [np.array(cC.to_rgb(i)) for i in c]  #list with RGB numpy array
                else:   #if4, else: the color cannot be RBG or RGBA
                    raise ValueError('Invalid list or array of rgb')
                #end if4
            else:  #if2, else: if more dimention
                raise ValueError('The rgb or rgba values must be contained in a 1D or 2D list or array')
            #end if2
        #end try1
    #end if1

    return rgbs

def _is_number(s):
    """
    Check if *c* is a number (from
    http://stackoverflow.com/questions/354038/how-do-i-check-if-a-string-is-a-number-in-python)
    Parameters
    ----------
    c: variable
    output
    ------
    true if c is a number
    false otherwise
    """
    try:
        float(s) # for int, long and float
    except ValueError:
        return False
    return True

def _check_alpha(alpha, n):
    """
    Check if alpha has one or n elements and if they are numberics and between 0 and 1
    Parameters
    ----------
    alpha: number or list/tuple/numpy array of numbers
        values to check
    output
    ------
    alpha: list of numbers 
        if all elements numberics and between 0 and 1
    """
    alpha = np.array(alpha).flatten()  #convert alpha to a flattened array
    if(alpha.size == 1):  #if1: alpha is one element
        if(_is_number(alpha) == False or alpha < 0 or alpha > 1):
            raise ValueError("'alpha' must be a float with value between 0 and 1, included") 
        else:
            alpha = [alpha for i in range(n)]  #replicate the alphas len(colors) times
    elif(alpha.size==n):  #if1, else: if alpha is composed of len(colors) elements
        try:  #check if all alphas are numbers
            alpha+1 
        except TypeError:
            raise ValueError("All elements of alpha must be a float with value between 0 and 1, included") 
        else:
            if((alpha < 0).any() or (alpha > 1).any()):
                raise ValueError("'alpha' must be a float with value between 0 and 1, included") 
    else:  #if1, else: if none of the previous cases
        raise ValueError("Alpha must have either one element or as many as 'colors'")
    #end if1
    return alpha

def colorAlpha_to_rgb(colors, alpha, bg='w'):
    """
    Given a Matplotlib color and a value of alpha, it returns 
    a RGB color which mimic the RGBA colors on the given background
    Parameters
    ----------
    colors: Matplotlib color (documentation from matplotlib.colors.colorConverter.to_rgb), 
        list/tuple/numpy array of colors
        Can be an *RGB* or *RGBA* sequence or a string in any of
        several forms:
        1) a letter from the set 'rgbcmykw'
        2) a hex color string, like '#00FFFF'
        3) a standard name, like 'aqua'
        4) a float, like '0.4', indicating gray on a 0-1 scale
        if *color* is *RGBA*, the *A* will simply be discarded.
    alpha: float [0,1] or list/tuple/numpy array with len(colors) elements
        Value of alpha to mimic. 
    bg: Matplotlib color (optional, default='w')
        Color of the background. Can be of any type shown in *color*
    output
    ------
    rgb: *RGB* color 
    example
    -------
    import mimic_alpha as ma
    print(ma.colorAlpha_to_rgb('r', 0.5))
    >>> [array([ 1. ,  0.5,  0.5])]
    print(ma.colorAlpha_to_rgb(['r', 'g'], 0.5)) 
    >>> [array([ 1. ,  0.5,  0.5]), array([ 0.5 ,  0.75,  0.5 ])]
    print(ma.colorAlpha_to_rgb(['r', 'g'], [0.5, 0.3])) 
    >>> [array([ 1. ,  0.5,  0.5]), array([ 0.7 ,  0.85,  0.7 ])]
    print(ma.colorAlpha_to_rgb(['r', [1,0,0]], 0.5)) 
    >>> [array([ 1. ,  0.5,  0.5]), array([ 1. ,  0.5,  0.5])]
    print( ma.colorAlpha_to_rgb([[0,1,1], [1,0,0]], 0.5) ) 
    >>> [array([ 0.5,  1. ,  1. ]), array([ 1. ,  0.5,  0.5])]
    print(ma.colorAlpha_to_rgb(np.array([[0,1,1], [1,0,0]]), 0.5)) 
    >>> [array([ 0.5,  1. ,  1. ]), array([ 1. ,  0.5,  0.5])]
    print(ma.colorAlpha_to_rgb(np.array([[0,1,1], [1,0,0]]), 0.5, bg='0.5')) 
    >>> [array([ 0.25,  0.75,  0.75]), array([ 0.75,  0.25,  0.25])]
    """

    colors = _to_rgb(colors)  #convert the color and save in a list of np arrays
    bg = np.array(cC.to_rgb(bg))   #convert the background

    #check if alpha has 1 or len(colors) elements and return a list of len(color) alpha 
    alpha = _check_alpha(alpha, len(colors))  
    #interpolate between background and color 
    rgb = [(1.-a) * bg + a*c for c,a in zip(colors, alpha)]

    return rgb


def cmap(cmap_name, alpha, bg="w", set_under=None, set_over=None,
         set_bad=None, out_cmap_name=None):
    """
    Generate an RGB colormap from a given mpl cmap and alpha value.
    Parameters
    ----------
    cmap_name: String
       A standard Matplotlib colormap name:
       http://matplotlib.org/examples/color/colormaps_reference.html
    alpha: Float
       Value of alpha to mimic in range [0,1].
    bg: Matplotlib color
       Color of the background.
    out_cmap_name: String
       Name of the returned colormap.
    set_under: Matplotlib color
       Set color to be used for low out-of-range values.
    set_over: Matplotlib color
       Set color to be used for high out-of-range values.
    set_bad: Matplotlib color
       Set color to be used for masked values.
    Output
    ------
    ma_cmap: :class:`matplotlib.colors.Colormap`
       A colormap instance that mimics an RGBA standard cmap.
    Notes
    -----
    This code is based on the make_cmap() program written by Chris Slocum:
      http://schubert.atmos.colostate.edu/~cslocum/custom_cmap.html
    Example
    -------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import mimic_alpha as ma
    >>> plt.ion()
    >>> # Make a gradient image:
    >>> gradient = np.linspace(0, 1, 50)
    >>> image = np.repeat(np.atleast_2d(gradient), repeats=2, axis=0)
    >>> # Compare contourf() plots without alpha, with alpha, and mimic-alpha:
    >>> plt.figure(0)
    >>> plt.clf()
    >>> plt.subplots_adjust(0.1, 0.1, 0.95, 0.95, hspace=0.4)
    >>> ax = plt.subplot(411)
    >>> ax.set_title("Standard 'hot' colormap")
    >>> cs = plt.contourf(image, levels=gradient, cmap="hot")
    >>> ax.set_xticklabels([""])
    >>> ax = plt.subplot(412)
    >>> mahot = ma.cmap("hot", 1.0)
    >>> ax.set_title("Mimic-alpha 'hot' colormap with alpha=1.0")
    >>> cs = plt.contourf(image, levels=gradient, cmap=mahot)
    >>> ax.set_xticklabels([""])
    >>> ax = plt.subplot(413)
    >>> ax.set_title("Standard 'hot' colormap with alpha=0.5")
    >>> cs = plt.contourf(image, levels=gradient, cmap="hot", alpha=0.5)
    >>> ax.set_xticklabels([""])
    >>> ax = plt.subplot(414)
    >>> mahot = ma.cmap("hot", 0.5)
    >>> ax.set_title("Mimic-alpha 'hot' colormap with alpha=0.5")
    >>> cs = plt.contourf(image, levels=gradient, cmap=mahot)
    >>> # Compare outputs when saved as a postscript file:
    >>> plt.savefig("mimic_alpha_hot.ps")
    """
    # Read input cmap:
    input_cmap = plt.cm.get_cmap(cmap_name)
    ncolors = input_cmap.N

    position = np.linspace(0, 1, ncolors)
    # Convert RGBA colors from cmap into RGB:
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos in position:
        r, g, b = colorAlpha_to_rgb(input_cmap(pos), alpha, bg)[0]
        cdict['red'  ].append((pos, r, r))
        cdict['green'].append((pos, g, g))
        cdict['blue' ].append((pos, b, b))

    # Set output colormap name:
    if out_cmap_name is None:
        out_cmap_name = cmap_name + "_{0:.1f}".format(alpha)
    # mimic-alpha colormap:
    ma_cmap = mplc.LinearSegmentedColormap(out_cmap_name, cdict, 256)

    # Set mimic-alpha colors for masked and out-of-range values:
    if set_under is not None:
        RGBunder = colorAlpha_to_rgb(set_under, alpha, bg)[0]
        ma_cmap.set_under(RGBunder)
    if set_over is not None:
        RGBover  = colorAlpha_to_rgb(set_over,  alpha, bg)[0]
        ma_cmap.set_over(RGBover)
    if set_bad is not None:
        RGBbad   = colorAlpha_to_rgb(set_bad,   alpha, bg)[0]
        ma_cmap.set_bad(RGBbad)

    return ma_cmap