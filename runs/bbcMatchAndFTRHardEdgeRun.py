import numpy as np
import matplotlib.pyplot as plt
from momentSolver.MomentSolver import MomentSolver,MomentSolverUtility
from momentSolver.PlottingUtility import PlottingUtility
# Magnet parameters & Opt parameters
from systems.bbcMatchingAndFTRHardEdge.magnetParameters import lattice

msu = MomentSolverUtility()
pu = PlottingUtility()

# initial values
betax,betay = 0.629, 0.0629
alphax,alphay = 0,0
emitx,emity = 53e-6,5.3e-6
#initCond = msu.GetInitialConditions(betax,betay,alphax,alphay,emitx,emity)

xprimesq=2.118983100854736858e-05
xprimeyprime=0
xsq=8.307505759310390044e-06
xxprime=0
xy=0
xyprime=0
yprimesq=2.104241921605963170e-05
ysq=8.320693194604818684e-08

yxprime=0

yyprime=-2.091128042711822252e-09
#initCond = msu.GetInitialConditions(betax,betay,alphax,alphay,emitx,emity)
#initCond=msu.GetInitialConditionsFromWarp( xsq, xprimesq, ysq, yprimesq, xy, xxprime, yyprime, xyprime, yxprime, xprimeyprime)

initCond=msu.GetInitialConditionsFromWarp( xsq, xprimesq, ysq, yprimesq, xy, xxprime, yyprime, xyprime, yxprime, xprimeyprime)#
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
mom.Run(verbose=True)

# plot results
pu.PlotEnv(mom)
plt.show()