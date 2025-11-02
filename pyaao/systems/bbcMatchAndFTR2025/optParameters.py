from MomentSolver import OptimizationUtility
import numpy as np

ou = OptimizationUtility()

# run optimization on just FTR quad strength and nothing else in the lattice
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
    paramarray = [1.00689455,
                1.00738497,
                0.955682,
                0.96906425,
                0.98753207,
                1.11503013,
                0.98753207,
                0.9915573]
    
# run an optimization on all the magnets
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

    # 15 optimization parameters. Initialize everything to 1.0 to start
    paramarray = np.ones(15)     

params = ou.getParamObj(paramarray, parammapping)