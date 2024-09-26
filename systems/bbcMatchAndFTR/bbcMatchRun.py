import numpy as np
import matplotlib.pyplot as plt
from momentSolver import MomentSolver,MomentSolverUtility
from momentSolver import PlottingUtility

# import your lattice
from magnetParameters import lattice

# initialize helpful functions
msu = MomentSolverUtility()
pu = PlottingUtility()

# initial values
betax,betay = 0.629, 0.0629
alphax,alphay = 0,0
emitx,emity = 53e-6,5.3e-6
slitx, slity = 0.01, 0.001 # slit size in meters
# betax,betay = (slitx)**2/(3.0*emitx), (slity)**2/(3.0*emity)
initCond = msu.GetInitialConditions(betax,betay,alphax,alphay,emitx,emity)

xprimesq = 2.188e-5
xprimeyprime = 3.498e-8
xsq = 8.30e-6
xxprime = -4.45e-8
xy = -1.96e-10
xyprime = 1.49e-8
yprimesq = 2.10e-5
ysq = 8.32e-8
yxprime = -1.59e-9
yyprime = -2.09e-9

# physics settings
energy = 5e3 # eV
current = 2.9e-3 # Amps
pipeRadius = 0.0 # meters

# sim parameters
zInterval = (0, 2.0) # meters
stepSize = 0.0001 # step size

# initialize moment solver
mom = MomentSolver(
    lattice, 
    initialConditions=initCond, 
    energy=energy, 
    current=current, 
    pipeRadius=pipeRadius, 
    zInterval=zInterval, 
    stepSize=stepSize
)

# run the moment equations over the given distance
mom.Run(verbose=True)

# plot results
pu.PlotEnv(mom)
plt.show()

# if you just want the data to do your own analysis/plotting:
Ydata = mom.y
zpos = mom.z

# zpos is a list of z points, e.g. 0 meters - > 0.33 meters
# Ydata is an 10 by N array containing all the moment values, Q+,Q-,Qx,P+,P-,Px,E+,E-,Ex,L as a function of z position

# e.g. if I wanted to plot Q+ vs z for my run, I can do something like:
# plt.plot(zpos, Ydata[0,:])

# take a look at the PlotEnv function being called above to see how we make that plot etc...