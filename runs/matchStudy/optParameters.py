import numpy as np
from momentSolver import OptimizationUtility

ou = OptimizationUtility()

if True:

    paramarray = np.array([ 
        1.0,
        1.0,
        1.0,
        1.0,
    ])  

    # order to map param list to param object
    parammapping = [['dbdx'],['dbdx',],['dbdx'],['dbdx']] 

    quad1 = {}
    quad1['dbdx'] = paramarray[0]

    quad2 = {}     
    quad2['dbdx'] = paramarray[1]

    quad3 = {}    
    quad3['dbdx'] = paramarray[2]

    quad4 = {}    
    quad4['dbdx'] = paramarray[3]    

params = ou.getParamObj(paramarray, parammapping)