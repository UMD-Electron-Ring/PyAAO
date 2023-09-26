# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 12:12:07 2023

@author: lpocher
"""

import numpy as np
import matplotlib.pyplot as plt
from momentSolver import MomentSolver,OptimizationUtility,MomentSolverUtility,PlottingUtility
                
# Magnet parameters & Opt parameters
from magnetParametersHardEdge import lattice
from optParameters import params,parammapping

# restrictions on the optimization. This gets applied during each gradient step
def setRestrictions(momObj, paramArray):

    # right now false, so there is no restriction. If we want to add one just make this true and put in the logic
    # if False:
        # We want quad 1 & 3 to have same strengths, so in our case, have the same optimization values        
        # Thus, in our optimization parameters, quad 1 strength and quad 3 strength were parameters 0 and 2
        # paramArray[0] = paramArray[2]# 7/21 I want the 1st and the 3rd to be the same

    # only make the solenoid able to be moved back not forward for LSE
    # if(paramArray[3]<1.0):
    #    paramArray[3] = 1.0

    # other examples of restrictions you can do:
    # If you wanted to restrict the position of magnets
    # say you want quad 1 and quad 2 zstart positions to be atleast 0.18meters away, we can do some logic like:
    # if ( paramArray[3]['zstart']*mom.latticeDefault[3]['zstart'] - mom.lattice[2]['zend'] < 0.05164 ):
    #     paramArray[3]['zstart'] = (0.05164 + mom.lattice[2]['zstart'])/mom.latticeDefault[3]['zstart']
    # Logic is a bit more complicated here since the optimization parameters are scalars, we need to get the actual z positions and multiply by scalars
    # Let me know if you want to do stuff like this and I can help
    # I am leaving the below as kinda more examples
    # if False:
    #     length = 0.12 # 12 cm
    #     quad3Center = momObj.lattice[2]['zstart'] * paramArray[6] + momObj.lattice[2]['length'] / 2.0
    #     paramArray[9] = (length + quad3Center) / (momObj.latticeDefault[-1]['zstart'])
    #     #paramArray[9] = 0.8
    #     #print( "POS: " + str(quad3Center) + "|" + str(paramArray[9] * momObj.lattice[3]['zstart']) )        

    return paramArray, momObj

msu = MomentSolverUtility()
ou = OptimizationUtility(setRestrictions) # restrictions function handle gets set here
pu = PlottingUtility()

# initial values
betax,betay = 0.629, 0.0629
alphax,alphay = 0,0
emitx,emity = 53e-6,5.3e-6
initCond = msu.GetInitialConditions(betax,betay,alphax,alphay,emitx,emity)

# physics settings
energy = 5e3 # [eV]
current = 0.0e-3# [Amps]
pipeRadius = 0.0 # [meters] , for image charges effect on pipe walls, zero ignores the effect

# sim parameters
# zInterval = (0.0, 1.80114) 
# zInterval = (0.0, 0.5) # got an interestign results for 0.5 m run length
zInterval = (0.0,.6489) # mechanical end of the matching section from slit
stepSize = .1e-3

print('Setting up initial lattice')
mom = MomentSolver(lattice, energy=energy, current=current, pipeRadius=pipeRadius, zInterval=zInterval, stepSize=stepSize)
mom.initialMoments = initCond

# lattice gets initialize with the initial magnet values from the magnetParameters file.
# if we want to then update these with initial optimization parameters, in case the opt parameters aren't all 1.0 starting, we do it like this:
mom.UpdateLattice(params = params)
pu.printLattice(mom) # nice tabular printout 

# run moment equations and plot whatever we have
mom.Run()
pu.PlotEnv(mom, title='Initial solution')

###############################################################

# run moment solver. Ctrl-c to interrupt
mom, an_h, gamma_h, f_h, fp_h, df_h  = ou.runOptimization(mom, params, parammapping)

# save results
# an_h is an history of the value of optimization parameters over every iteration step
# gamma_h is the history of the gamma parameters every iteration step
# f_h is the history of the FoM every iteration step
# fp_h is the history of each component of the FoM every iteration step, this is used mostly for debugging to see how much each components contributes to optimization
# df_h is gradient used in the gradient descent optimization, every iteration step
# I usually save data like this incase I need to reload old optimization values/parameters and start a new run from those values
#np.save('run1.npy',[an_h,gamma_h,f_h,fp_h,df_h],dtype=object)

# make your own plots here if you want
if True:
    plt.figure()
    plt.plot(mom.z, mom.y[0,:]) # plots Q+ vs z
    x2 = mom.y[0,:] + mom.y[1,:] # Q+ + Q- = <x^2>
    plt.plot(mom.z, x2) # plots <x^2> vs z

# print final results 
print(an_h[-1]) # print the last best set of optimization parameters
an_h[-1],mom = setRestrictions(mom, an_h[-1]) # apply the restriction function to these last set of opt parameters
mom.UpdateLattice( params = ou.getParamObj(an_h[-1],parammapping) ) # update the lattice with the last set of opt parameters
pu.printLattice(mom) # print out the lattice values

# plot stuff again, post optimization
mom.Run()
pu.PlotEnv(mom, title='Optimized solution')

plt.show()
