import numpy as np
from momentSolver import OptimizationUtility

ou = OptimizationUtility()

if True:

    paramarray = np.array([ 
        0.97349958,  
        1.17730812,  
        1.22577242,  
        1.6338717,   
        1.13054699,  
        1.31214369,
        1.51125034,  
        1.09008688,  
        1.01156225,  
        0.8,         
        -0.94121212
    ])    

    # order to map param list to param object
    parammapping = [['zstart','rotation','dbdx'],['zstart','rotation','dbdx'],['zstart','rotation','dbdx'],['zstart','dbdx']]    

    params = ou.getParamObj(paramarray, parammapping)
