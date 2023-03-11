from momentSolver.MomentSolver import OptimizationUtility
import numpy as np

ou = OptimizationUtility()

# run with all parameters being used
if False:
    paramarray = np.ones(11)  

    # order to map param list to param object
    parammapping = [['zstart','rotation','dbdx'],['zstart','rotation','dbdx'],['zstart','rotation','dbdx'],['zstart','dbdx']]    

    quad1 = {}
    quad1['zstart'] = paramarray[0]
    quad1['rotation'] = paramarray[1]
    quad1['dbdx'] = paramarray[2]

    quad2 = {}
    quad2['zstart'] = paramarray[3]
    quad2['rotation'] = paramarray[4]
    quad2['dbdx'] = paramarray[5]

    quad3 = {}
    quad3['zstart'] = paramarray[6]
    quad3['rotation'] = paramarray[7]
    quad3['dbdx'] = paramarray[8]

    sol = {}
    sol['zstart'] = paramarray[9]
    sol['dbdx'] = paramarray[10]    

# run with only magnet strength, no movement or rotations
if True:
    paramarray = np.ones(4) 

    paramarray = np.array([
        1.00920768, 
        1.01666039, 
        1.03008797, 
        0.9948411, 
        1.0
        ])

    # order to map param list to param object
    #parammapping = [['rotation','dbdx'],['rotation','dbdx'],['rotation','dbdx'],['dbdx']]
    parammapping = [['dbdx'],['dbdx','zstart'],['dbdx'],['dbdx', 'zstart']]

    quad1 = {}
    quad1['dbdx'] = paramarray[0]

    quad2 = {}
    quad2['dbdx'] = paramarray[1]
    quad2['zstart'] = paramarray[2]    

    quad3 = {}
    quad3['dbdx'] = paramarray[3]

    sol = {}
    sol['dbdx'] = paramarray[4]
    sol['zstart'] = paramarray[5]

# run with magnet strength and rotations, no movement except solenoid
if False:
    paramarray = np.array([
        1.06700521, 
        1.19518164, 
        0.99439876, 
        1.1122049,  
        1.01345775, 
        1.04320169,
        0.95774436
    ])        

    #[1.03295908 1.04032733 1.02197807 1.04688933 1.01706542 1.03425943 0.99670053]    
    #[1.12394992 1.27467979 1.03634494 1.1696602  1.05025807 1.10642427 0.94499498]

    # order to map param list to param object
    #parammapping = [['rotation','dbdx'],['rotation','dbdx'],['rotation','dbdx'],['dbdx']]
    parammapping = [['rotation','dbdx'],['rotation','dbdx'],['rotation','dbdx'],['dbdx']] 

params = ou.getParamObj(paramarray, parammapping)