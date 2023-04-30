import random
import numpy as np
import matplotlib.pyplot as plt
from momentSolver.MomentSolver import MomentSolver,OptimizationUtility,MomentSolverUtility
from momentSolver.MomentSolver import OptimizationUtility
from momentSolver.PlottingUtility import PlottingUtility
# Magnet parameters & Opt parameters
from systems.bbcMatchingAndFTRHardEdge.magnetParameters import lattice
#from systems.bbcMatchingAndFTRHardEdge.optParameters import params,parammapping

msu = MomentSolverUtility()
ou = OptimizationUtility()
pu = PlottingUtility()

## Taken from OptParameters

ou = OptimizationUtility()

# run optimization all the following values:
# strengths of all 4 matching quads
# positions, rotations, and strength of FTR triplet
# position and strength of solenoid
if False:
    # order of paramters in lattice
    parammapping = [
        ['dbdx'], # match quad 1 strength
        ['dbdx'], # match quad 2 strength
        ['dbdx'], # match quad 3 strength
        ['dbdx'], # match quad 4 strength
        ['zstart','rotation','dbdx'], # FTR quad 1
        ['zstart','rotation','dbdx'], # FTR quad 2
        ['zstart','rotation','dbdx'], # FTR quad 3
        ['zstart','dbdx'] # solenoid
    ]      

    paramarray = np.ones(15) 

# run optimization on just FTR quad strength
if True:
    # order of paramters in lattice
    parammapping = [
        ['None'], # match quad 1 
        ['None'], # match quad 2
        ['None'], # match quad 3 
        ['None'], # match quad 4
        ['dbdx'], # FTR quad 1
        ['dbdx'], # FTR quad 2
        ['dbdx'], # FTR quad 3
        ['None'] # solenoid
    ]      

    # we have 8 things in our mapping ,so we need 8 initial values for them (1)
    # Even though the 'None' categories will be ignored, we need to just initialize some value for them
    paramarray = np.ones(8)

# run optimization on just FTR quad strength and solenoid strength + position
if False:
    # order of paramters in lattice
    parammapping = [
        ['None'], # match quad 1 
        ['None'], # match quad 2
        ['None'], # match quad 3 
        ['None'], # match quad 4
        ['dbdx'], # FTR quad 1
        ['dbdx'], # FTR quad 2
        ['dbdx'], # FTR quad 3
        ['zstart', 'dbdx'] # solenoid
    ]      

    paramarray = np.ones(10)    

params = ou.getParamObj(paramarray, parammapping)

random.seed(223)
for i,value in enumerate(paramarray):
   start=.98*value# array notationa not needed for modifying
   stop=1.02*value
   paramarray[i]=np.round(np.random.uniform(start, stop), 6)
for x in paramarray:
   print("Random number within 2% of " + str(x) + " = ")
   start = .98 * x# need brackets otherwise function
   stop = 1.02 * x
   output = np.round(np.random.uniform(start, stop,10), 6)
   print(str(output) + "\n")

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

# FoM to use
useFoMForMatching = False

# initial values
betax,betay = 0.629, 0.0629
alphax,alphay = 0,0
emitx,emity = 53e-6,5.3e-6

xprimesq=2.118983100854736858e-05
xprimeyprime=3.498629953673319272e-08
xsq=8.307505759310390044e-06
xxprime=-4.456914897925242239e-08
xy=-1.967466871186281683e-10
xyprime=1.497600674
yprimesq=2.104241921605963170e-05
ysq=8.320693194604818684e-08

yxprime=-1.590008975507985115e-09

yyprime=-2.091128042711822252e-09

initCond = msu.GetInitialConditions(betax,betay,alphax,alphay,emitx,emity)
#initCond = msu.GetInitialConditionsWarp(betax,betay,alphax,alphay,emitx,emity)

# physics settings
energy = 5e3 # [eV]
current = 2.9e3 # [Amps]
pipeRadius = 0.0 # [meters] , for image charges effect on pipe walls, zero ignores the effect

# sim parameters
zInterval = (0, 2.5) # 1.422 is about where the solenoid starts
stepSize = 0.0001




print('Setting up initial lattice')
mom = MomentSolver(lattice, energy=energy, current=current, pipeRadius=pipeRadius, zInterval=zInterval, stepSize=stepSize)
mom.initialMoments = initCond
mom.UpdateLattice(params = params)
pu.printLattice(mom)

# plot initial stuff
mom.Run()
pu.PlotEnv(mom, title='Initial solution')
mom.zInterval = (0, 2.5)
mom.Run()
pu.PlotEnv(mom, title='Initial solution over longer distance solenoid')
mom.zInterval = zInterval

# run moments and adjoint equations
print('Running Mom. Eqn.')
for b in range(10):
# run moment solver. Ctrl-c to interrupt
    mom, an_h, gamma_h, f_h, fp_h, df_h  = ou.runOptimization(mom, params, parammapping)

# main loop 

    #Plot of the figure of merit versus iternation number in the optmization loop (roughly)
    plt.figure()
    plt.plot(range(len(f_h)),f_h)# supposed to be plotting figure of merit versus iteration number in the optimization lop, is this the same as the x xis in the standard plots generated by this code
    #plt.show()
    #Plot of Error Plots
    plt.figure()
    xr=mom.y[0,:]+mom.y[1,:]
    yr=mom.y[0,:]-mom.y[1,:]
    Ratio=np.divide(xr,yr)
    Difference=np.substract(xr,yr)
    Error=np.sqrt(np.divide(Difference**2,xr))
    plt.plot(mom.z, Ratio, label='$<Ratio>$')
    #plt.show()
    plt.plot(mom.z, Error, label='$<Error>$')
    plt.show()
#print final results
    an_h[-1],mom = setRestrictions(mom, an_h[-1])
    mom.UpdateLattice( params = ou.getParamObj(an_h[-1],parammapping) )
    pu.printLattice(mom)

    # plot stuff again
    pu.PlotEnv(mom, title='Optimized solution')
    mom.zInterval = (0, 2.5)
    mom.Run()
    pu.PlotEnv(mom, title='Optimized solution over longer distance solenoid')
    mom.zInterval = zInterval

    plt.show()
    #print(b,"after plot commands")