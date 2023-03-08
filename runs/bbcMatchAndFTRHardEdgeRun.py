from momentSolver.MomentSolver import MomentSolver
import numpy as np
import matplotlib.pyplot as plt

###################################
# Magnet parameters
from systems.bbcMatchingAndFTRHardEdge.magnetParameters import lattice
###################################

# initial beam conditions
def GetInitialConditions(betax, betay, alphax, alphay, emitx, emity):
    xrms = np.sqrt(betax * emitx)
    yrms = np.sqrt(betay * emity)
    Q_plus = 0.5*(xrms**2 + yrms**2)
    Q_minus = 0.5*(xrms**2 - yrms**2)
    Q_x = 0.0
    P_plus = (-alphax - alphay)
    P_minus = (-alphax + alphay)
    P_x = 0
    E_plus = 0.5 * ( 2*(emitx**2 + 0.25 * (P_plus + P_minus)**2)/(Q_plus + Q_minus) + 2*(emity**2 + 0.25 * (P_plus - P_minus)**2)/(Q_plus - Q_minus) )
    E_minus = 0.5 * ( 2*(emitx**2 + 0.25 * (P_plus + P_minus)**2)/(Q_plus + Q_minus) - 2*(emity**2 + 0.25 * (P_plus - P_minus)**2)/(Q_plus - Q_minus) )
    E_x = 0
    L = 0
    phi = 0
    return np.array([Q_plus,Q_minus,Q_x,P_plus,P_minus,P_x,E_plus,E_minus,E_x,L,phi])

# initial values
betax,betay = 0.629, 0.0629
alphax,alphay = 0,0
emitx,emity = 53e-6,5.3e-6
initCond = GetInitialConditions(betax,betay,alphax,alphay,emitx,emity)

# physics settings
energy = 5e3 # eV
current = 0.0e-3 # Amps
pipeRadius = 0.0 # meters

# sim parameters
zInterval = (0, 2.5) # meters
stepSize = 0.0001 # step size

mom = MomentSolver(
    lattice, 
    initialConditions=initCond, 
    energy=energy, 
    current=current, 
    pipeRadius=pipeRadius, 
    zInterval=zInterval, 
    stepSize=stepSize
)

# run the moment equations
z, y, ksol, kquad = mom.Run(verbose=True)

xr = y[0,:] + y[1,:] #  <x^2> = Q+ + Q-
yr = y[0,:] - y[1,:] #  <y^2> = Q+ - Q-
plt.figure()
plt.plot(z,np.sqrt(xr), label='$<x>$')
plt.plot(z,np.sqrt(yr), label='$<y>$')
plt.plot(z, ksol * 6e-5, c='m') # plot scaled magnet plots in the same plot so we know where the magnets are in the lattice
plt.plot(z,kquad * 3e-6, c='k') # plot scaled magnet plots in the same plot so we know where the magnets are in the lattice
plt.xlabel('Z [m]')
plt.ylabel('Beam moment [m]')
plt.grid(True)
plt.legend()

plt.show()