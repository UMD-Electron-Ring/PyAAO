import numpy as np
from momentSolver import OptimizationUtility

ou = OptimizationUtility()

if True:

    paramarray = np.array([ 
        1.0,
        1.0,
        1.0,
    ])  

    # order to map param list to param object
    parammapping = [['dbdx'],['dbdx',],['dbdx']] 

params = ou.getParamObj(paramarray, parammapping)