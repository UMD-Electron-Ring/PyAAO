import numpy as np
import matplotlib.pyplot as plt
from momentSolver.MomentSolver import MomentSolver,MomentSolverUtility, OptimizationUtility
from momentSolver.PlottingUtility import PlottingUtility
# Magnet parameters & Opt parameters
from systems.bbcMatchingAndFTRHardEdge.magnetParameters import lattice

msu = MomentSolverUtility()
pu = PlottingUtility()
ou = OptimizationUtility()

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
current = 2.9e-3 # Amps
pipeRadius = 0.0 # meters

# sim parameters
zInterval = (0, 2.5) # meters
stepSize = 0.0001 # step size

####### Setup Monte Carlo Runs
# order of paramters in lattice
parammapping = [
    ['dbdx'], # match quad 1 
    ['dbdx'], # match quad 2
    ['dbdx'], # match quad 3 
    ['dbdx'], # match quad 4
    ['dbdx'], # FTR quad 1
    ['dbdx'], # FTR quad 2
    ['dbdx'], # FTR quad 3
    ['dbdx'] # solenoid
]      
paramarray = np.array([
    1.18712312, 
    1.09061345, 
    0.94645688, 
    1.00492277, 
    0.94393287, 
    1.14678234,
    0.94393287, 
    0.98379653
])
params = ou.getParamObj(paramarray, parammapping)

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
#plt.ylim((-1,10))

mom.UpdateLattice(params = params)

# run the moment equations
mom.Run(verbose=True)
# plot results
pu.PlotEnv(mom)
#plt.ylim((-1,10))

plt.show()