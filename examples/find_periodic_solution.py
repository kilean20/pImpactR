# there already is a much faster matching code named 3D ENV 
# this example is just an illustration of parallel optimization with IMPACT 

import os
import numpy as np
import cPickle as pickle
import shutil 
import pIMPACT as pm

#%%
beam, lattice = pm.readIMPACT('FODO.in')

beam['energy']=300E6
beam['n_particles']=500
beam['current']=0.0
#if space-charge considered uncomment followings
#beam['n_particles']=3000
#beam['current']=0.04


QuadIndex=[]
QuadStrength=[]
for i in range(len(lattice)):
    if lattice[i]['type']=='quad':
        QuadIndex.append(i)
        QuadStrength.append(lattice[i]['B1'])
        lattice[i]['n_sckick']=int(np.ceil(lattice[i]['length']*40))
    if lattice[i]['type']=='drift':  
        lattice[i]['n_sckick']=int(np.ceil(lattice[i]['length']*1))
        #if space-charge considered uncomment followign
        #lattice[i]['n_sckick']=int(np.ceil(lattice[i]['length']*7))

#%%
def objFunc(arg): 
    for i in QuadIndex:
        lattice[i]['B1']=(-1)**(i/2)*arg[0]
    pm.twiss2beam(beam, arg[1], -arg[2], 0.26E-6, 
                         arg[1], arg[2], 0.26E-6, 
                         30.0, 0.0, 0.314)
    
    target = pm.opt.id_generator()  # generage random directory name
    while os.path.exists(target):  
        target = pm.opt.id_generator()
    shutil.copytree('origin', target) # copy working directory to random directroy
    # In this example there is no data file needed in './origin' 
    
    os.chdir(target) # cd to the randome directory and
    
    pm.writeIMPACT('test.in',beam,lattice)
    pm.run() # run impact there
    
    X=pm.readBeamSize('x')
    Y=pm.readBeamSize('y')
    betx, alfx, enx, bety, alfy, eny, betz, alfz, enz = pm.readOpticsAtEnd()
    
    # optimize average rms size to 1.5 mm
    obj1 = np.sum( (np.array(X)*1E3-1.5)**2 + (np.array(Y)*1E3-1.5)**2 )
    # periodic condition
    obj2 = (betx - arg[1])**4 + (5.0*alfx + 5.0*arg[2])**4 +\
           (bety - arg[1])**4 + (5.0*alfy - 5.0*arg[2])**4
    os.chdir('..')
    shutil.rmtree(target)
    return obj1 + 10*obj2
#%% run optim
bounds = [(7.0,14.0), (6.0,12.0), (0.9,2.1)]
result=pm.opt.differential_evolution(objFunc, bounds, ncore=24, popsize=8, 
                                     disp=True, polish=False, maxtime=60*1) 
                                     # stop running at maximum 1 min
# save current population of optimization 
with open('result.data','wb') as fp:
    pickle.dump(result,fp)
  
#%% resume optimization until converge. 
# Ability to resume optimization is very useful especially for NERSC debug mode
while True:
    with open('result.data','rb') as fp:
        previous_result = pickle.load(fp)
    result = pm.opt.differential_evolution(objFunc, bounds, ncore=24, 
                                           prev_result=previous_result, 
                                           disp=True, polish=False, maxtime=60*1)
    if hasattr(result,'x'): 
        break                            


#%%
# print optimization result and save in directory ./print_result
def print_result(arg): 
    for i in QuadIndex:
        lattice[i]['B1']=(-1)**(i/2)*arg[0]
    pm.twiss2beam(beam, arg[1], -arg[2], 0.26E-6, 
                         arg[1], arg[2], 0.26E-6, 
                         30.0, 0.0, 0.314)
    
    target = 'print_result'
    shutil.copytree('origin', target)
    os.chdir(target)      
    
    pm.writeIMPACT('test.in',beam,lattice)
    pm.run() # run impact there
    
    X=pm.readBeamSize('x')
    Y=pm.readBeamSize('y')
    betx, alfx, enx, bety, alfy, eny, betz, alfz, enz = pm.readOpticsAtEnd()
    pm.plot.rms()
    os.chdir('..')

    return betx,alfx,enx,bety,alfy,eny,betz,alfz,enz, max(X),max(Y)

print print_result(result.x)
