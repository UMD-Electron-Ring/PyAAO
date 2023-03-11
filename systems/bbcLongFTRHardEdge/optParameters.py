from momentSolver.MomentSolver import OptimizationUtility
import numpy as np

ou = OptimizationUtility()

# run with 
if True:
    paramarray = np.ones(4)  

    quad1 = {}
    quad1['dbdx'] = paramarray[0]

    quad2 = {}
    quad2['dbdx'] = paramarray[1]

    quad3 = {}
    quad3['dbdx'] = paramarray[2]

    sol = {}
    sol['dbdx'] = paramarray[3]    

# order to map param list to param object
parammapping = [['dbdx'],['dbdx'],['dbdx'],['dbdx']]
params = ou.getParamObj(paramarray, parammapping)