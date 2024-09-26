# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 11:36:42 2023

@author: lpocher
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os # for directory creation etc.
from momentSolver import MomentSolver,OptimizationUtility,MomentSolverUtility,PlottingUtility
                
# Magnet parameters & Opt parameters
from magnetParametersHardEdge import lattice
from optParameters import params,parammapping

plt.close('all')

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

# plotting stuff
ls = ["solid","dashed","dotted","dashdot"]
mark = ['+', '.', 'o', '*']
pnames = ["Init_Sol",
          "Opt_Sol",
          "Q_Sol",
          "P_Sol",
          "E_Sol",
          "Params",
          "FoM",
          "FoM_i"
          ]

# get current working directort
cwd = os.getcwd()
dname = cwd + "\3FTR"

# units for easy conversion
um = 1.0e-6; mm = 1.0e-3; us = um; ms = mm

# restrictions on simulation/optimization
mSep = 0.12 # minimum separation between quadroles
q1mDist = 0.155 # minimum distance quad 1 center can get to slit  0.121539 / 0.172339
quadmount = 21.25*0.0254 # quad mount length
leff = 0.05164 # effective length of quadrupoles in meters for hardedge
lqact = 1.16*2.0*0.0254 # actual physical length of quadrupole
sollen = 1.37 # physical length of the solenoid
zSol = q1mDist - lqact/2.0 + quadmount # starting location of physical solenoid
zL = zSol + sollen/2.0 # total length of simulation halfway through sol.

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
    mSep = 0.12 # minimum separation between quadroles
    q1mDist = 0.155 # minimum distance quad 1 can get to slit  0.121539 / 0.172339
    quadmount = 21.25*0.0254 # quad mount length
    leff = 0.05164 # effective length of quadrupoles in meters
    lqact = 1.16*2.0*0.0254 # actual physical length of quadrupole
    ls = 2.0 # length of solenoid field is 2 meters
    sollen = 1.37
    zSol = q1mDist - lqact/2.0 + quadmount # starting location of physical solenoid
    zL = zSol + 1.37/2.0 # total length of simulation halfway through sol.
    
    q1loc = paramArray[0]*mom.latticeDefault[0]['zstart']
    q2loc = paramArray[2]*mom.latticeDefault[1]['zstart']
    q3loc = paramArray[4]*mom.latticeDefault[2]['zstart']
    # q4loc = paramArray[6]*mom.latticeDefault[3]['zstart']
    q1cen = q1loc + leff/2.0
    q2cen = q2loc + leff/2.0
    q3cen = q3loc + leff/2.0
    # q4cen = q4loc + leff/2.0
    
    # solloc = paramArray[6]*mom.latticeDefault[3]['zstart']
    # solpstart = solloc # physical start of the long sol.
    
    # minimum distance
    if ( q1cen < q1mDist ):
        paramArray[0] = (q1mDist-leff/2.0)/(mom.latticeDefault[0]['zstart'])
        
    if ( q2cen < (mSep + q1loc) ):
        paramArray[2] = (mSep + q1loc)/(mom.latticeDefault[1]['zstart'])
        
    if ( q3cen < (mSep + q2loc) ):
        paramArray[4] = (mSep + q2loc)/(mom.latticeDefault[2]['zstart'])
        
    # if ( q4cen < (mSep + q3loc) ):
    #     paramArray[6] = (mSep + q3loc)/(mom.latticeDefault[3]['zstart'])

    # maximum distance
    if ( q1cen > (zSol - 2.0*mSep - leff/2.0 ) ):
        paramArray[0] = (zSol - 3*mSep - leff)/(mom.latticeDefault[0]['zstart'])
        
    if ( q2cen > (zSol - 1.0*mSep - leff/2.0 ) ):
        paramArray[2] = (zSol - 2.0*mSep - leff)/(mom.latticeDefault[1]['zstart'])
        
    if ( q3cen > (zSol - 0.0*mSep - leff/2.0 ) ):
        paramArray[4] = (zSol - 1.0*mSep - leff)/(mom.latticeDefault[2]['zstart'])
        
    # if ( q4cen > (zSol - 0.0*mSep - leff/2.0 ) ):
    #     paramArray[6] = (zSol - 0.0*mSep - leff)/(mom.latticeDefault[3]['zstart'])
       
    # if ( solpstart < zSol):
    #     paramArray[6] = (zSol)/(mom.latticeDefault[3]['zstart'])

    return paramArray, momObj


# for loop iteration on the things to iterate on
Niter = 11
alpha0, alphaf = -1.0, 1.0 # rate of (con)divergence of the beam
zFinal0, zFinalf = zSol, 0.99# optimization locations in meters
alphaVec = np.linspace(alpha0,alphaf,Niter)
zFinalVec = np.linspace(zFinal0,zFinalf,Niter)

# get str to save files to/make directory for saving things
dname = "FTR45_3QuadHardedge/"
if (os.path.exists(dname)==0):
    os.mkdir(dname)
    

for alpha in alphaVec:
    for zF in zFinalVec:
        
        plt.close('all') # close all plots
        print("\nAlphaX = %.2f \t zF = %.2f\n"%(alpha,zF))
        
        # initial values
        emitx,emity = 20.0e-6, 2.0e-6 # what's the emittance after the slit?
        slitx, slity = 0.01, 0.001 # slit size in x
        betax,betay = slitx**2/(3.0*emitx), slity**2/(3.0*emity)
        alphax = alpha
        alphay = alphax*(betay/betax)
       
        
        # physics settings
        energy = 5e3 # [eV]
        current = 2.6e-3# [Amps]
        pipeRadius = 0.0 # [meters] , for image charges effect on pipe walls, zero ignores the effect
        
        # sim parameters
        zInterval = (0.0, zF) # mechanical end of the matching section from slit
        stepSize = .1e-3
        
        # file name str
        sdname = "emitx%.2e_betax%.2f_Ek%.2e_Ib%.2e_alphax%.2f_zF%.2f/"%(emitx,betax,energy,current,alphax,zF)
        pname = r"$\varepsilon=$%.2e $\beta_x=$%.2f $E_k=$%.2e $I_b=$%.2e $\alpha_x=$%.2f $z_F=$%.2f"%(emitx,betax,energy,current,alphax,zF)
        if (os.path.exists(dname + sdname) == 0):
            os.mkdir(dname + sdname)
        if (os.path.exists(dname + sdname + "Parameters.npy") == 1):
            continue
            
        # make objects and get IC
        msu = MomentSolverUtility()
        ou = OptimizationUtility(setRestrictions) # restrictions function handle gets set here
        pu = PlottingUtility()
        initCond = msu.GetInitialConditions(betax,betay,alphax,alphay,emitx,emity)
        
        print('Setting up initial lattice')
        mom = MomentSolver(lattice, energy=energy, current=current, pipeRadius=pipeRadius, zInterval=zInterval, stepSize=stepSize)
        mom.initialMoments = initCond
        
        # lattice gets initialize with the initial magnet values from the magnetParameters file.
        # if we want to then update these with initial optimization parameters, in case the opt parameters aren't all 1.0 starting, we do it like this:
        mom.UpdateLattice(params = params)
        pu.printLattice(mom) # nice tabular printout 
        
        # run moment equations and plot whatever we have
        mom.Run()
        pu.PlotEnv(mom, title=r"Init. Sol." + pname)
        
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
        
        # print final results 
        print(an_h[-1]) # print the last best set of optimization parameters
        an_h[-1],mom = setRestrictions(mom, an_h[-1]) # apply the restriction function to these last set of opt parameters
        mom.UpdateLattice( params = ou.getParamObj(an_h[-1],parammapping) ) # update the lattice with the last set of opt parameters
        pu.printLattice(mom) # print out the lattice values
        
        # zIntervalEnd = (0.0, an_h[-1][6]*mom.latticeDefault[3]['zstart'] + sollen)
        zIntervalEnd = (0.0, zF + sollen)
        # make a second moment object with the lattice parameters of the optimized lattice
        momTotalLen = MomentSolver(lattice, energy=energy, current=current, pipeRadius=pipeRadius, zInterval=zIntervalEnd, stepSize=stepSize)
        momTotalLen.initialMoments = initCond
        
        # lattice gets initialize with the initial magnet values from the magnetParameters file.
        # if we want to then update these with initial optimization parameters, in case the opt parameters aren't all 1.0 starting, we do it like this:
        momTotalLen.UpdateLattice(params = ou.getParamObj(an_h[-1],parammapping) )
        pu.printLattice(momTotalLen) # nice tabular printout 
        
        # plot stuff again, post optimization
        momTotalLen.Run()
        pu.PlotEnv(momTotalLen, title=r"Opt. Sol." + pname)
        
        plt.figure() # plot Q's
        plt.plot(momTotalLen.z, momTotalLen.y[0,:]/um,linestyle='-',marker=mark[0],label=r'$Q_+$')
        plt.plot(momTotalLen.z, momTotalLen.y[1,:]/um,linestyle='--',marker=mark[1],label=r'$Q_-$')
        plt.plot(momTotalLen.z, momTotalLen.y[2,:]/um,linestyle='-.',marker=mark[2],label=r'$Q_x$')
        plt.grid(True)
        plt.xlabel(r"$z$ [m]")
        plt.ylabel(r"$Q$ [mm$^2$]")
        plt.legend(loc="best")
        plt.tight_layout()
        
        plt.figure() # plot P's
        plt.plot(momTotalLen.z, momTotalLen.y[3,:]/um,linestyle='-',marker=mark[0],label=r'$P_+$')
        plt.plot(momTotalLen.z, momTotalLen.y[4,:]/um,linestyle='--',marker=mark[1],label=r'$P_-$')
        plt.plot(momTotalLen.z, momTotalLen.y[5,:]/um,linestyle='-.',marker=mark[2],label=r'$P_x$')
        plt.plot(momTotalLen.z, momTotalLen.y[9,:]/um,linestyle=':',marker=mark[2],label=r'$L$')
        plt.grid(True)
        plt.xlabel(r"$z$ [m]")
        plt.ylabel(r"$P$ and $L$ [mm mrad]")
        plt.legend(loc="best")
        plt.tight_layout()
        
        plt.figure()
        plt.plot(momTotalLen.z, momTotalLen.y[6,:]/um,linestyle='-',marker=mark[0],label=r'$E_+$')
        plt.plot(momTotalLen.z, momTotalLen.y[7,:]/um,linestyle='--',marker=mark[1],label=r'$E_-$')
        plt.plot(momTotalLen.z, momTotalLen.y[8,:]/um,linestyle='-.',marker=mark[2],label=r'$E_x$')
        plt.grid(True)
        plt.xlabel(r"$z$ [m]")
        plt.ylabel(r"$E$ [mrad$^2$]")
        plt.legend(loc="best")
        plt.tight_layout()
        
        # Plot History of Parameter Variation
        pu.plotParams(an_h)
        
        # plto Figure of Merit
        pu.plotfom(f_h,fp_h)
        
        # save figures (all 8 of them as .pngs)
        fig_nums = plt.get_fignums()
        figs = [plt.figure(n) for n in fig_nums]
        for i,fig in enumerate(figs):
            fig.savefig(dname + sdname + pnames[i] + ".png", dpi=200,transparent=True)
        
        # save data
        np.save(dname + sdname + "Parameters.npy",np.array([an_h,gamma_h,f_h,fp_h,df_h],dtype=object))
        
        plt.show() # show figs
        
        # delete variables so they can be resued
        del msu, ou, pu, mom, momTotalLen
        
# EOF