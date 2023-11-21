import numpy as np
from momentSolver import OptimizationUtility

ou = OptimizationUtility()

if True:

    paramarray = np.array([    
        1.0, # q1 loc
        1.0, # q1 str
        1.0, # q2 loc
        1.0, # q2 str
        1.0, # q3 loc
        1.0, # q3 str
        1.0, # q4 loc
        1.0, # q4 str
        1.0, # sol loc
        1.0  # sol str
    ])    
    
    
    # paramarray = np.array([    
    #     1.0, # q1 loc
    #     1.0, # q1 str
    #     1.0, # q2 loc
    #     1.0, # q2 str
    #     1.0, # q3 loc
    #     1.0, # q3 str
    # ])  

    
    # paramarray = np.array([    
    #     1.0, # q1 loc
    #     1.0, # q1 str
    #     1.0, # q2 loc
    #     1.0, # q2 str
    #     1.0, # q3 loc
    #     1.0, # q3 str
    #     1.0, # q4 loc
    #     1.0  # q4 str
    # ])  
    
    # paramarray = np.array([    
    #     1.0, # q1 loc
    #     1.0, # q1 str
    #     1.0, # q2 loc
    #     1.0, # q2 str
    #     1.0, # q3 loc
    #     1.0, # q3 str
    #     1.0, # q4 loc
    #     1.0, # q4 str
    #     1.0  # sol str
    # ])  

    # order to map param list to param object
    # parammapping = [['zstart','dbdx'],['zstart','dbdx'],['zstart','dbdx']]
    # parammapping = [['zstart','dbdx'],['zstart','dbdx'],['zstart','dbdx'],['zstart','dbdx']]
    parammapping = [['zstart','dbdx'],['zstart','dbdx'],['zstart','dbdx'],['zstart','dbdx'],['zstart','dbdx']]

    params = ou.getParamObj(paramarray, parammapping)
