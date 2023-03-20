import numpy as np
import matplotlib.pyplot as plt
from momentSolver.MomentSolver import MomentSolver,MomentSolverUtility
from momentSolver.PlottingUtility import PlottingUtility

# Magnet parameters
from systems.bbcLongFTRHardEdge.magnetParameters import lattice

msu = MomentSolverUtility()
pu = PlottingUtility()

# initial values
betax,betay = 0.629, 0.0629
alphax,alphay = 0,0
emitx,emity = 53e-6,5.3e-6
initCond = msu.GetInitialConditions(betax,betay,alphax,alphay,emitx,emity)

# physics settings
energy = 5e3 # eV # 5000eV or 5keV or 0.005MeV
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
mom.Run(verbose=True)

# plot results
pu.PlotEnv(mom)
plt.show()