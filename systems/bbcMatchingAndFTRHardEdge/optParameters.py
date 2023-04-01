from momentSolver.MomentSolver import OptimizationUtility
import numpy as np

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