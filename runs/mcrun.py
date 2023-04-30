import numpy as np
import matplotlib.pyplot as plt
from momentSolver.MomentSolver import MomentSolver,OptimizationUtility,MomentSolverUtility
from momentSolver.PlottingUtility import PlottingUtility
# Magnet parameters & Opt parameters
from systems.bbcMatchingAndFTRHardEdge.magnetParameters import lattice

# restrictions on the optimization. This gets applied during each gradient step
def setRestrictions(momObj, paramArray):

    # FTR quad 1 & 3 symmetric strength restriction
    # We want quad 1 & 3 to have same strengths, so in our case, have the same optimization values
    if True:
        # Figure out what order the dbdx strength parameters were in our paramArray
        # If we are running using the second set of parameters in the optParameters.py file, then quad1 & 3 strength is 
        # inside index number 4 & 6:
        # Let us just make quad 1 dbdx the same as quad 3 dbdx
        paramArray[4] = paramArray[6]

    return paramArray, momObj

msu = MomentSolverUtility()
ou = OptimizationUtility(setRestrictions)
pu = PlottingUtility()

# initial values
betax,betay = 0.629, 0.0629
alphax,alphay = 0,0
emitx,emity = 53e-6,5.3e-6

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
initCond=msu.GetInitialConditionsFromWarp( xsq, xprimesq, ysq, yprimesq, xy, xxprime, yyprime, xyprime, yxprime, xprimeyprime)#constant

# physics settings
energy = 5e3 # [eV]
current = 2.9e-3 # [Amps]
pipeRadius = 0.0 # [meters] , for image charges effect on pipe walls, zero ignores the effect

# sim parameters
zInterval = (0, 1.35) # 1.422 is about where the solenoid starts
stepSize = 0.0001

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
    1.000, 
    1.000, 
    1.000, 
    1.000, 
    1.000, 
    1.000,
    1.000, 
    1.000,
])

# generate a list of initial paramarray values to use.

# say 10 different values
paramarrays = [] # list of a set of paramarray
Nparameters = 8
Nruns = 1
import random
for i in range(Nruns):
    # can do something better here, e.g. uniform start etc..
    random.seed(223)
    placeholder=[]# defining
    for i,value in enumerate(paramarray):
        start=.98*value# array notationa not needed for modifying
        stop=1.02*value
        placeholder.append(np.round(np.random.uniform(start, stop), 6))
   
   
    #randValues = (np.random.rand(Nparameters) * 0.2 + 0.9) # will give you a random number between 0.9 -> 1.10 , so +/- 10%
    #paramarrays.append( paramarray * randValues )
    paramarrays.append( placeholder )# lists of different parameter sets
# for i in range(Nruns):
#     # can do something better here, e.g. uniform start etc..
#     randValues = (np.random.rand(Nparameters) * 0.2 + 0.9) # will give you a random number between 0.9 -> 1.10 , so +/- 10%
#     paramarrays.append( paramarray * randValues )

# params = ou.getParamObj(paramarray, parammapping)
#######

###############################################################

# arrays to say stuff in
# e.g. an_h_mc stands for an_history_montecarlo . Each run will return a an_h, so with X monte carlos we will have X an_h that we will store in an_h_mc
an_h_mc = []
gamma_h_mc = []
f_h_mc = []
fp_H_mc = []
df_h_mc = []
mom_mc = []

# run optimization for each set of paramarrays in our lsit
for ii,paramarrayChoice in enumerate(paramarrays):

    params = ou.getParamObj(paramarrayChoice, parammapping)
    print('Setting up initial lattice')
    mom = MomentSolver(lattice, energy=energy, current=current, pipeRadius=pipeRadius, zInterval=zInterval, stepSize=stepSize)
    mom.initialMoments = initCond
    mom.UpdateLattice(params = params)
    pu.printLattice(mom)

    # initial plot if you want
    if True:
        # plot initial stuff
        mom.Run()
        pu.PlotEnv(mom, title='Initial solution, run #' + str(ii))   

    # run moment solver. Ctrl-c to interrupt
    mom, an_h, gamma_h, f_h, fp_h, df_h  = ou.runOptimization(mom, params, parammapping, maxSteps=500)#originally 1000

    # print final results
    print(an_h[-1])
    an_h[-1],mom = setRestrictions(mom, an_h[-1])
    mom.UpdateLattice( params = ou.getParamObj(an_h[-1],parammapping) )
    pu.printLattice(mom)

    # save final results
    # each mc run is a separate file
    np.save('runs/run'+str(ii)+'.npy',[mom.z,mom.y,an_h,gamma_h,f_h,fp_h,df_h])

    # grab all the data from the run and store it.
    an_h_mc.append( an_h )
    gamma_h_mc.append( gamma_h )
    f_h_mc.append( f_h )
    fp_H_mc.append( fp_h )
    df_h_mc.append( df_h )
    mom_mc.append( mom )  

    if True:
        # plot stuff at the end
        mom.Run()
        pu.PlotEnv(mom, title='Optimized solution, run #' + str(ii))

# save final results
# each mc run is a separate file
np.save('runs/runTotal.npy',[an_h_mc,gamma_h_mc,f_h_mc,fp_H_mc,df_h_mc])        


plt.show()