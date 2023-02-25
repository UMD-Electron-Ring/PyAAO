from momentSolver.MomentSolver import MomentSolver
import numpy as np
import matplotlib.pyplot as plt

###################################
# Magnet parameters
from systems.bbcInitialMatch.magnetParameters import lattice
###################################

###################################
# Opt parameters
from systems.bbcInitialMatch.optParameters import params, initialConditions
###################################

# physics settings
energy = 5e3 # eV
current = 0.0e-3 # Amps
pipeRadius = 0.0 # meters

# sim parameters
zInterval = (0, 0.5) # meters
stepSize = 0.0001 # step size

mom = MomentSolver(lattice, energy=energy, current=current, pipeRadius=pipeRadius, zInterval=zInterval, stepSize=stepSize)
mom.initialMoments = initialConditions
mom.UpdateLattice(params = params)

z, y, ksol, kquad = mom.Run(verbose=True)
#zadj, yadj, _, _ = mom.RunAdjoint(verbose=True)

xr = y[0,:] + y[1,:] #  <x^2> = Q+ + Q-
yr = y[0,:] - y[1,:] #  <y^2> = Q+ - Q-
plt.figure()
plt.plot(z,xr, label='beam size X')
plt.plot(z,yr, label='beam size Y')
plt.plot(z,ksol * 1e-6) # plot scaled magnet plots in the same plot so we know where the magnets are in the lattice
plt.plot(z,kquad * 1e-8) # plot scaled magnet plots in the same plot so we know where the magnets are in the lattice
plt.xlabel('Z [m]')
plt.ylabel('Beam moment [m^2]')
plt.grid(True)
plt.legend()

plt.show()