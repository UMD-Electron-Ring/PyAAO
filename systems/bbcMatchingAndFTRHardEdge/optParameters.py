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

if True:
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

    paramarray = np.ones(15) 

params = getParamObj(paramarray, parammapping)