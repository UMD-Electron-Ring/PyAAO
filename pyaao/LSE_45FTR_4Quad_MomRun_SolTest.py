# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 15:09:06 2023

@author: lpocher
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 10:40:39 2023

@author: lpocher
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from momentSolver import MomentSolver,OptimizationUtility,MomentSolverUtility,PlottingUtility
                
# Magnet parameters & Opt parameters
from magnetParametersHardEdgeSol import lattice
from optParametersSol import params,parammapping

pFlag = False # plotting flag for lots of diagnostic plots

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

    return paramArray, momObj

msu = MomentSolverUtility()
ou = OptimizationUtility(setRestrictions) # restrictions function handle gets set here
pu = PlottingUtility()

# physics settings
energy = 5e3 # [eV]
current = 2.6e-3# [Amps]
# current = 0.0e-3# [Amps] # great solution at about z = 0.55 for no space charge!
pipeRadius = 0.0 # [meters] , for image charges effect on pipe walls, zero ignores the effect

# sim parameters
# zInterval = (0.0, 1.80114) 
# zInterval = (0.0, 0.5) # got an interestign results for 0.5 m run length
zF = zSol + sollen
zInterval = (0.0, zF) # mechanical end of the matching section from slit
stepSize = .1e-3

# initial values
emitx,emity = 20.0e-6, 2.0e-6 # what's the emittance after the slit?
slitx, slity = 0.01, 0.001 # slit size in x
betax,betay = (slitx)**2/(3.0*emitx), (slity)**2/(3.0*emity)
alphax = 0.0
alphay = alphax*(betay/betax)
initCond = msu.GetInitialConditions(betax,betay,alphax,alphay,emitx,emity)

mom = MomentSolver(lattice, energy=energy, current=current, pipeRadius=pipeRadius, zInterval=zInterval, stepSize=stepSize)
mom.initialMoments = initCond

pname = r"$\varepsilon=$%.2e $\beta_x=$%.2f $E_k=$%.2e $I_b=$%.2e $\alpha_x=$%.2f $z_F=$%.2f"%(emitx,betax,energy,current,alphax,zF)

# lattice gets initialize with the initial magnet values from the magnetParameters file.
# if we want to then update these with initial optimization parameters, in case the opt parameters aren't all 1.0 starting, we do it like this:
mom.UpdateLattice(params = params)
pu.printLattice(mom) # nice tabular printout 

# run moment equations and plot whatever we have
mom.Run()
pu.PlotEnv(mom, r"Mom. Sol. for $B_{LS}=%.2f$ G "%(mom.latticeDefault[3]['dbdx'].GetValue(zSol + 0.01)*10000) + pname )

if pFlag:
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