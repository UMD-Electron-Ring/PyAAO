# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 12:12:07 2023

@author: lpocher
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from momentSolver import MomentSolver,OptimizationUtility,MomentSolverUtility,PlottingUtility
                
# Magnet parameters & Opt parameters
from magnetParametersHardEdge import lattice
from optParameters import params,parammapping

# for plotting
csfont = {"fontname":"Times New Roman"}
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["text.usetex"] = True
# colors for plotting things
colors =  plt.rcParams["axes.prop_cycle"].by_key()["color"]
# plotting stuff
title_size = 28
label_size = 24
leg_size   = 16
tick_size = 16
GR = 1.618
COLOR = "black"
mpl.rcParams["text.color"] = COLOR
mpl.rcParams["axes.labelcolor"] = COLOR
mpl.rcParams["xtick.color"] = COLOR
mpl.rcParams["ytick.color"] = COLOR

ls = ["solid","dashed","dotted","dashdot"]
mark = ['+', '.', 'o', '*']

um = 1.0e-6; mm = 1.0e-3; us = um; ms = mm

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
    
    #### enforce minimin and maximum distances between magnets
    # get locations and set minimum separation between the magnets
    q1loc = paramArray[0]*mom.latticeDefault[0]['zstart']
    q2loc = paramArray[2]*mom.latticeDefault[1]['zstart']
    q3loc = paramArray[4]*mom.latticeDefault[2]['zstart']
    mSep = 0.12 # minimum separation
    q1mDist = 0.01 # minimum distance from the slit/start of simulatiom
    zL = .6489 # final extent of simulation
    leff = 0.05164 # effective length of quadrupoles in meters
    
    # minimum distance
    if ( q1loc < q1mDist ):
        paramArray[0] = q1mDist/(mom.latticeDefault[0]['zstart'])
        
    if ( q2loc < (mSep + q1loc) ):
        paramArray[2] = (mSep + q1loc)/(mom.latticeDefault[1]['zstart'])
        
    if ( q3loc < (q2loc + mSep) ):
        paramArray[4] = (mSep + q2loc)/(mom.latticeDefault[2]['zstart'])

    # maximum distance
    if ( q1loc > (zL - 2*mSep - leff ) ):
        paramArray[0] = (zL - 2*mSep - leff )/(mom.latticeDefault[0]['zstart'])
        
    if ( q2loc > (zL - mSep - leff ) ):
        paramArray[2] = (zL - mSep - leff )/(mom.latticeDefault[1]['zstart'])
        
    if ( q3loc > (zL - leff ) ):
        paramArray[4] = (zL - leff )/(mom.latticeDefault[2]['zstart'])

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
emitx,emity = 20.0e-6, 2.0e-6 # what's the emittance after the slit?
slitx, slity = 0.01, 0.001 # slit size in x
betax,betay = slitx**2/(3.0*emitx), slity**2/(3.0*emity)
alphax = 0.0
alphay = alphax*(betay/betax)
initCond = msu.GetInitialConditions(betax,betay,alphax,alphay,emitx,emity)

# physics settings
energy = 5e3 # [eV]
current = 2.6e-3# [Amps]
pipeRadius = 0.0 # [meters] , for image charges effect on pipe walls, zero ignores the effect

# sim parameters
# zInterval = (0.0, 1.80114) 
# zInterval = (0.0, 0.5) # got an interestign results for 0.5 m run length
zF = 0.7
zInterval = (0.0,zF) # mechanical end of the matching section from slit
stepSize = .1e-3

# plot names 
pname = r"$\varepsilon=$%.2e $\beta_x=$%.2f $E_k=$%.2e $I_b=$%.2e $\alpha_x=$%.2f $z_F=$%.2f"%(emitx,betax,energy,current,alphax,zF)

print('Setting up initial lattice')
mom = MomentSolver(lattice, energy=energy, current=current, pipeRadius=pipeRadius, zInterval=zInterval, stepSize=stepSize)
mom.initialMoments = initCond

# lattice gets initialize with the initial magnet values from the magnetParameters file.
# if we want to then update these with initial optimization parameters, in case the opt parameters aren't all 1.0 starting, we do it like this:
mom.UpdateLattice(params = params)
pu.printLattice(mom) # nice tabular printout 

# run moment equations and plot whatever we have
mom.Run()
pu.PlotEnv(mom, title='Init. Sol. ' + pname)

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
    plt.figure() # plot Q's
    plt.plot(mom.z, mom.y[0,:]/um,linestyle='-',marker=mark[0],label=r'$Q_+$')
    plt.plot(mom.z, mom.y[1,:]/um,linestyle='--',marker=mark[1],label=r'$Q_-$')
    plt.plot(mom.z, mom.y[2,:]/um,linestyle='-.',marker=mark[2],label=r'$Q_x$')
    plt.grid(True)
    plt.xlabel(r"$z$ [m]")
    plt.ylabel(r"$Q$ [mm$^2$]")
    plt.legend(loc="best")
    plt.tight_layout()
    
    plt.figure() # plot P's
    plt.plot(mom.z, mom.y[3,:]/um,linestyle='-',marker=mark[0],label=r'$P_+$')
    plt.plot(mom.z, mom.y[4,:]/um,linestyle='--',marker=mark[1],label=r'$P_-$')
    plt.plot(mom.z, mom.y[5,:]/um,linestyle='-.',marker=mark[2],label=r'$P_x$')
    plt.plot(mom.z, mom.y[9,:]/um,linestyle=':',marker=mark[2],label=r'$L$')
    plt.grid(True)
    plt.xlabel(r"$z$ [m]")
    plt.ylabel(r"$P$ and $L$ [mm mrad]")
    plt.legend(loc="best")
    plt.tight_layout()
    
    plt.figure()
    plt.plot(mom.z, mom.y[6,:]/um,linestyle='-',marker=mark[0],label=r'$E_+$')
    plt.plot(mom.z, mom.y[7,:]/um,linestyle='--',marker=mark[1],label=r'$E_-$')
    plt.plot(mom.z, mom.y[8,:]/um,linestyle='-.',marker=mark[2],label=r'$E_x$')
    plt.grid(True)
    plt.xlabel(r"$z$ [m]")
    plt.ylabel(r"$E$ [mrad$^2$]")
    plt.legend(loc="best")
    plt.tight_layout()

# print final results 
print(an_h[-1]) # print the last best set of optimization parameters
an_h[-1],mom = setRestrictions(mom, an_h[-1]) # apply the restriction function to these last set of opt parameters
mom.UpdateLattice( params = ou.getParamObj(an_h[-1],parammapping) ) # update the lattice with the last set of opt parameters
pu.printLattice(mom) # print out the lattice values

# plot stuff again, post optimization
mom.Run()
pu.PlotEnv(mom, title='Opt. Sol. ' + pname)

# Plot History of Parameter Variation
pu.plotParams(an_h)

# plto Figure of Merit
pu.plotfom(f_h,fp_h)

plt.show()
