import numpy as np
import matplotlib.pyplot as plt
from momentSolver import MomentSolver,MomentSolverUtility
from momentSolver import PlottingUtility

# import your lattice
#from magnetParameters import lattice
from systems.bbcMatchAndFTR2025.magnetParameters import lattice

# initialize helpful functions
msu = MomentSolverUtility()
pu = PlottingUtility()

liam_data_path = '/mnt/c/Users/levon/Downloads/Moment_Table_07-10-2025.txt'
liam_data = np.loadtxt(liam_data_path, delimiter=',', skiprows=1)

# initial values
# betax,betay = 0.629, 0.0629
# alphax,alphay = 0,0
# emitx,emity = 53e-6,5.3e-6
# initCond = msu.GetInitialConditions(betax,betay,alphax,alphay,emitx,emity)

# copy initial conditions straight from Liam data
initCond = np.array([
    (liam_data[0,1] + liam_data[0,4]) * 0.5, # Q+
    (liam_data[0,1] - liam_data[0,4]) * 0.5, # Q-
    liam_data[0,7], # Qx
    liam_data[0,3] + liam_data[0,6], # P+
    liam_data[0,3] - liam_data[0,6], # P-
    liam_data[0,10] + liam_data[0,8], # Px
    liam_data[0,2] + liam_data[0,5], # E+
    liam_data[0,2] - liam_data[0,5], # E-
    2.0 * liam_data[0,9], # Ex
    liam_data[0,8] - liam_data[0,10], # L
    0.0 # phi
])

# physics settings
energy = 5e3 # eV
current = 0.491e-3 # Amps
pipeRadius = 0.0 # meters

# sim parameters
zInterval = (0, 1.5) # meters
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
plt.plot(liam_data[:,0], np.sqrt(liam_data[:,1])*1e3,label='liam_x',color='k',linestyle='--')
plt.plot(liam_data[:,0], np.sqrt(liam_data[:,4])*1e3,label='liam_y',color='k',linestyle='--')
plt.legend()
plt.show()

# if you just want the data to do your own analysis/plotting:
Ydata = mom.y
zpos = mom.z

# zpos is a list of z points, e.g. 0 meters - > 0.33 meters
# Ydata is an 10 by N array containing all the moment values, Q+,Q-,Qx,P+,P-,Px,E+,E-,Ex,L as a function of z position

# e.g. if I wanted to plot Q+ vs z for my run, I can do something like:
# plt.plot(zpos, Ydata[0,:])

# take a look at the PlotEnv function being called above to see how we make that plot etc...