import numpy as np
import matplotlib.pyplot as plt
from momentSolver.MomentSolver import MomentSolver,OptimizationUtility,MomentSolverUtility
from momentSolver.PlottingUtility import PlottingUtility
# Magnet parameters & Opt parameters
from systems.bbcMatchingAndFTRHardEdge.magnetParameters import lattice

msu = MomentSolverUtility()
ou = OptimizationUtility()
pu = PlottingUtility()

# parameter values in lattice
parammapping = [
    ['dbdx'], # Matching quad 1
    ['dbdx'], # Matching quad 2
    ['dbdx'], # Matching quad 3    
    ['dbdx'], # Matching quad 4      
    ['dbdx'], # FTR quad 1
    ['dbdx'], # FTR quad 2
    ['dbdx'], # FTR quad 3
    ['None'] # solenoid
] 
paramarray = np.ones(8) # set all starting opt values to 1.0
params = ou.getParamObj(paramarray, parammapping)

# initial values
betax,betay = 0.629, 0.0629
alphax,alphay = 0,0
emitx,emity = 53e-6,5.3e-6
initCond = msu.GetInitialConditions(betax,betay,alphax,alphay,emitx,emity)

# physics settings
energy = 5e3 # eV
current = 0.0e-3 # Amps
pipeRadius = 0.0 # meters

# sim parameters
zInterval = (0, 1.4) # meters
stepSize = 0.0001 # step size

# we want the original , before any perturbations are made
momOriginal = MomentSolver(
    lattice, 
    initialConditions=initCond, 
    energy=energy, 
    current=current, 
    pipeRadius=pipeRadius, 
    zInterval=zInterval, 
    stepSize=stepSize
)

# copy we will use to perturb opt parameters
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
mom.Run()
mom.RunAdjoint(verbose=True)
momOriginal.Run()
momOriginal.RunAdjoint(verbose=True)

##################
# get initial, pre perturbation, values
O_nopert,N_nopert = momOriginal.Omat, momOriginal.Nmat
fom_nopert,_,_ = momOriginal.GetFoM_And_DFoM() # FoM

# lets perturb each magnet by about 20 pts and in a range from 0.8 -> 1.20 x
Npts = 20
perturbPoints = np.linspace(0.80,1.20,Npts)
clrs = ['C0','C1','C2','C3','C4','C5','C6','C7'] # plotting colors

plt.figure()
# 7 quadrupole magnets,
magnetidx = [0,1,2,3,4,5,6] # loop through all 7 of our quadrupole strength opt parameters
for i in magnetidx:

    # define some arrays
    fom_val_true = []
    fom_val_linear = []
    xvals = []
    paramVals = np.copy(paramarray)

    print(str(i) + "/" + str(7))
    for ii in range(Npts):

        # perturb the parameter, update lattice
        paramVals[i] = paramarray[i] * perturbPoints[ii]
        mom.UpdateLattice( params = ou.getParamObj(paramVals,parammapping) )

        # calculate x axis k-strength value
        xvals.append( mom.lattice[i]['dbdx'].GetValue(0) - mom.latticeDefault[i]['dbdx'].GetValue(0) )

        # run env. eqn with new perturbed magnet
        mom.Run()
        
        # true FoM calculation
        mom.Run() # run equations with new pert lattice
        fomtmp,_,_ = mom.GetFoM_And_DFoM() # grab the new FoM
        fom_val_true.append( fom_nopert - fomtmp ) # [ FoM new - FoM original (unperturbed) ]

        # integral adjoint linear FoM calculation
        O_pert, N_pert = momOriginal.CalcON( mom.lattice ) # grab the O & N matrices as a result of the perturbed magnet
        O_pert2, N_pert2 = [], []
        # take a difference between the O & N matrix values before and after the perturbation
        for jj in range(len(O_pert)):
            O_pert2.append ( O_pert[jj] - O_nopert[jj] )
            N_pert2.append ( N_pert[jj] - N_nopert[jj] )
        fom_val_linear.append( momOriginal.CalcInt(O_pert2, N_pert2) ) # calculate the change in FoM via the adjoint integral

    # plot the results each in their own mini subplot
    plt.subplot(3,3,i+1)
    plt.plot(np.array(xvals), fom_val_true, color = clrs[i], linestyle='--')
    plt.plot(-np.array(xvals), fom_val_linear, color = clrs[i])

plt.show()