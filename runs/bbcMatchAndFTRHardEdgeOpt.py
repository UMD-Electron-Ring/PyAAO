import numpy as np
import matplotlib.pyplot as plt
from momentSolver.MomentSolver import MomentSolver,OptimizationUtility,MomentSolverUtility
from momentSolver.PlottingUtility import PlottingUtility
# Magnet parameters & Opt parameters
from systems.bbcInitialMatch.magnetParameters import lattice
from systems.bbcInitialMatch.optParameters import params,parammapping


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
current = 3.0e-3# [Amps]
pipeRadius = 0.0 # [meters] , for image charges effect on pipe walls, zero ignores the effect

# sim parameters
zInterval = (0, 0.6305) # 0.6305 is where I think the matching section ends? Double check with Santiago or tech note.
stepSize = 0.0001

print('Setting up initial lattice')
mom = MomentSolver(lattice, energy=energy, current=current, pipeRadius=pipeRadius, zInterval=zInterval, stepSize=stepSize)
mom.initialMoments = initCond
mom.UpdateLattice(params = params)
pu.printLattice(mom)

# plot initial stuff
mom.Run()
pu.PlotEnv(mom, title='Initial solution')

###############################################################

# run moment solver. Ctrl-c to interrupt
mom, an_h, gamma_h, f_h, fp_h, df_h  = ou.runOptimization(mom, params, parammapping, useMatchFoM=True)

# save results
#np.save('runs/run2.npy',[an_h,gamma_h,f_h,fp_h,df_h],dtype=object)

# print final results 
print(an_h[-1])
mom.UpdateLattice( params = ou.getParamObj(an_h[-1],parammapping) )
pu.printLattice(mom)

# plot stuff again
pu.PlotEnv(mom, title='Optimized solution')

plt.show()