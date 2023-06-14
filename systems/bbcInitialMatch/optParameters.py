from momentSolver.MomentSolver import OptimizationUtility
import numpy as np

ou = OptimizationUtility()

if True:
# run optimization on just FTR quad strength
    # order of paramters in lattice
    parammapping = [
        ['dbdx'], # match quad 1 
        ['dbdx'], # match quad 2
        ['dbdx'], # match quad 3 
        ['dbdx'] # match quad 4
    ]      

    # we have 8 things in our mapping ,so we need 8 initial values for them (1)
    # Even though the 'None' categories will be ignored, we need to just initialize some value for them
    paramarray = [ 1.0, 1.0, 1.0, 1.0]

params = ou.getParamObj(paramarray, parammapping)    