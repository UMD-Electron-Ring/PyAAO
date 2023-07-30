import numpy as np
import matplotlib.pyplot as plt
from momentSolver import MomentSolver,OptimizationUtility,MomentSolverUtility,PlottingUtility
                
# Magnet parameters & Opt parameters
from magnetParameters import lattice
from optParameters import params,parammapping

# restrictions on the optimization. This gets applied during each gradient step
def setRestrictions(momObj, paramArray):

    # FTR quad 1 & 3 symmetric strength restriction
    # We want quad 1 & 3 to have same strengths, so in our case, have the same optimization values
    if True:
        # Figure out what order the dbdx strength parameters were in our paramArray
        # If we are running using the second set of parameters in the optParameters.py file, then quad1 & 3 strength is 
        # inside index number 4 & 6:
        # Let us just make quad 1 dbdx the same as quad 3 dbdx
        paramArray[0] = paramArray[2]# 7/21 I want the 1st and the 3rd to be the same
        paramArray[1]=paramArray[0]*2 # I want the 2nd quad to be twice the 1st and the 3rd

    return paramArray, momObj

msu = MomentSolverUtility()
ou = OptimizationUtility()
pu = PlottingUtility()

# initial values
betax,betay = 0.629, 0.0629
alphax,alphay = 0,0
emitx,emity = 53e-6,5.3e-6
initCond = msu.GetInitialConditions(betax,betay,alphax,alphay,emitx,emity)

# physics settings
energy = 5e3 # [eV]
current = 0.0e-3# [Amps]
pipeRadius = 0.0 # [meters] , for image charges effect on pipe walls, zero ignores the effect

# sim parameters
zInterval = (0, 0.33) # 0.33 is where the 3 quads end
stepSize = 0.0001

print('Setting up initial lattice')
mom = MomentSolver(lattice, energy=energy, current=current, pipeRadius=pipeRadius, zInterval=zInterval, stepSize=stepSize)
mom.initialMoments = initCond
mom.UpdateLattice(params = params)
pu.printLattice(mom)

# plot initial stuff
mom.Run()
pu.PlotEnv(mom, title='Initial solution')
mom.zInterval = (0, 0.33)
mom.Run()
pu.PlotEnv(mom, title='Initial solution over longer distance solenoid')
mom.zInterval = zInterval

###############################################################

# run moment solver. Ctrl-c to interrupt
mom, an_h, gamma_h, f_h, fp_h, df_h  = ou.runOptimization(mom, params, parammapping)

# save results
#np.save('run1.npy',[an_h,gamma_h,f_h,fp_h,df_h],dtype=object)

# make your own plots here
plt.figure()
plt.plot(mom.z, mom.y[0,:]) # plots Q+ vs z
x2 = mom.y[0,:] + mom.y[1,:] # Q+ + Q- = <x^2>
plt.plot(mom.z, x2) # plots <x^2> vs z

# print final results 
print(an_h[-1])
an_h[-1],mom = setRestrictions(mom, an_h[-1])
mom.UpdateLattice( params = ou.getParamObj(an_h[-1],parammapping) )
pu.printLattice(mom)

# plot stuff again
pu.PlotEnv(mom, title='Optimized solution')
mom.zInterval = (0, 0.33)
mom.Run()
pu.PlotEnv(mom, title='Optimized solution over longer distance solenoid')
mom.zInterval = zInterval

plt.show()