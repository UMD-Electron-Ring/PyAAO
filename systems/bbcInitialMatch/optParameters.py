import numpy as np

def getParamArray(paramObj):
    arr = []
    for i,elem in enumerate(paramObj):
        for elemName in elem:
            arr.append(paramObj[i][elemName])
    return np.array(arr)

def getParamObj(paramArray):
    params = []
    numMags = len(paramArray)
    for i in range(numMags):
        tmp = {}
        tmp['dbdx'] = paramArray[i]
        params.append(tmp)
    return np.array(params)

if True:

    paramarray = np.array([ 
        1.0,
        1.0,
        1.0,
        1.0,
    ])    

    quad1 = {}
    quad1['dbdx'] = paramarray[0]

    quad2 = {}     
    quad2['dbdx'] = paramarray[1]

    quad3 = {}    
    quad3['dbdx'] = paramarray[2]

    quad4 = {}    
    quad4['dbdx'] = paramarray[3]    

params = np.array([ quad1, quad2, quad3, quad4])

# Initial conditions
Q_plus = 0.5*(2.2581**2*1e-6 + 0.2258**2*1e-6)
Q_minus = 0.0
Q_x = 0.0
P_plus = 0
P_minus = 0
P_x = 0
E_plus = (7.0855**2*1e-6 + 0.70855**2*1e-6)
E_minus = 0.0
E_x = 0
L = 0
phi = 0

initialConditions = np.array([Q_plus,Q_minus,Q_x,P_plus,P_minus,P_x,E_plus,E_minus,E_x,L,phi])