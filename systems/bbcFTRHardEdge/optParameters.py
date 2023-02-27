import numpy as np

def getParamArray(paramObj):
    arr = []
    for i,elem in enumerate(paramObj):
        for elemName in elem:
            arr.append(paramObj[i][elemName])
    return np.array(arr)

def getParamObj(paramArray, paramMapping):
    params = []
    cc = 0
    for ii,elem in enumerate(paramMapping):
        tmp = {}
        for jj,param in enumerate(elem):
            tmp[param] = paramArray[cc]
            cc += 1
        params.append(tmp)
    return np.array(params)   

# initial beam conditions
def GetInitialConditions(betax, betay, alphax, alphay, emitx, emity):
    xrms = np.sqrt(betax * emitx)
    yrms = np.sqrt(betay * emity)
    Q_plus = 0.5*(xrms**2 + yrms**2)
    Q_minus = 0.5*(xrms**2 - yrms**2)
    Q_x = 0.0
    P_plus = (-alphax - alphay)
    P_minus = (-alphax + alphay)
    P_x = 0
    E_plus = 0.5 * ( 2*(emitx**2 + 0.25 * (P_plus + P_minus)**2)/(Q_plus + Q_minus) + 2*(emity**2 + 0.25 * (P_plus - P_minus)**2)/(Q_plus - Q_minus) )
    E_minus = 0.5 * ( 2*(emitx**2 + 0.25 * (P_plus + P_minus)**2)/(Q_plus + Q_minus) - 2*(emity**2 + 0.25 * (P_plus - P_minus)**2)/(Q_plus - Q_minus) )
    E_x = 0
    L = 0
    phi = 0
    return np.array([Q_plus,Q_minus,Q_x,P_plus,P_minus,P_x,E_plus,E_minus,E_x,L,phi])

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
if False:
    paramarray = np.ones(4)  

    # order to map param list to param object
    #parammapping = [['rotation','dbdx'],['rotation','dbdx'],['rotation','dbdx'],['dbdx']]
    parammapping = [['dbdx'],['dbdx'],['dbdx'],['dbdx']]

    quad1 = {}
    #quad1['rotation'] = paramarray[0]
    quad1['dbdx'] = paramarray[0]

    quad2 = {}
    #quad2['rotation'] = paramarray[2]
    quad2['dbdx'] = paramarray[1]

    quad3 = {}
    #quad3['rotation'] = paramarray[4]
    quad3['dbdx'] = paramarray[2]

    sol = {}
    sol['dbdx'] = paramarray[3]

# run with magnet strength and rotations, no movement
if True:
    paramarray = np.array([
        1.0,
        1.02035703, 
        1.0,
        1.03109128, 
        1.0,
        1.02559988, 
        0.99755865
    ])

    # order to map param list to param object
    #parammapping = [['rotation','dbdx'],['rotation','dbdx'],['rotation','dbdx'],['dbdx']]
    parammapping = [['rotation','dbdx'],['rotation','dbdx'],['rotation','dbdx'],['dbdx']]

    quad1 = {}
    quad1['rotation'] = paramarray[0]
    quad1['dbdx'] = paramarray[1]

    quad2 = {}
    quad2['rotation'] = paramarray[2]
    quad2['dbdx'] = paramarray[3]

    quad3 = {}
    quad3['rotation'] = paramarray[4]
    quad3['dbdx'] = paramarray[5]

    sol = {}
    sol['dbdx'] = paramarray[6]    

if False:

    paramarray = np.ones(11)
    paramarray[-1] = -1.0

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

params = np.array([ quad1, quad2, quad3, sol])