from momentSolver.MomentSolver import MomentSolver
import numpy as np
import matplotlib.pyplot as plt

###################################
# Magnet parameters & Opt parameters
from systems.bbcFTRHardEdge.magnetParameters import lattice
from systems.bbcFTRHardEdge.optParameters import params,GetInitialConditions
###################################

# initial values
betax,betay = 0.698,0.698
alphax,alphay = 0,0
emitx,emity = 53e-6,5.3e-6
initCond = GetInitialConditions(betax,betay,alphax,alphay,emitx,emity)

# physics settings
energy = 5e3 # eV
current = 0.0e-3 # Amps
pipeRadius = 0.0 # meters

# sim parameters
zInterval = (0, 0.7) # meters
stepSize = 0.0001 # step size

mom = MomentSolver(lattice, energy=energy, current=current, pipeRadius=pipeRadius, zInterval=zInterval, stepSize=stepSize)
mom.UpdateLattice(params = params)
mom.initialMoments = initCond

z, y, ksol, kquad = mom.Run(verbose=True)
#zadj, yadj, _, _ = mom.RunAdjoint(verbose=True)

xr = y[0,:] + y[1,:] #  <x^2> = Q+ + Q-
yr = y[0,:] - y[1,:] #  <y^2> = Q+ - Q-
plt.figure()
plt.plot(z,xr, label='$<x^2>$')
plt.plot(z,yr, label='$<y^2>$')
plt.plot(z,-1 * ksol * 6e-7, c='m') # plot scaled magnet plots in the same plot so we know where the magnets are in the lattice
plt.plot(z,kquad * 3e-8, c='k') # plot scaled magnet plots in the same plot so we know where the magnets are in the lattice
plt.xlabel('Z [m]')
plt.ylabel('Beam moment [m^2]')
plt.grid(True)
plt.legend()

plt.show()