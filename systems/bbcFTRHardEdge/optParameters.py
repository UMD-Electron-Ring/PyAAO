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

    paramarray = np.array([
        1.02086247, 
        1.03050807, 
        1.02404281, 
        0.97240616,
        ]) 

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
if False:
    paramarray = np.array([
        1.03980134, 
        1.10362798, 
        0.99803625, 
        1.06436912, 
        1.00889097, 
        1.02682633,
        0.97279436])

    #[1.03295908 1.04032733 1.02197807 1.04688933 1.01706542 1.03425943 0.99670053]    

    # 3 mA
    #[1.08644173 1.25569886 0.99908952 1.14673349 1.01993958 1.05985367 0.95058752]    

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

# run with magnet strength and rotations, no movement except solenoid
if True:
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

if False:

    paramarray = np.array([
        1.00067974, 
        1.02192697, 
        1.01120623, 
        1.06126934, 
        1.01383525, 
        1.04485096,
        1.02811584, 
        1.00306504, 
        1.00029057,
        0.9774864,  
        0.99896648
    ])
    
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

params = getParamObj(paramarray, parammapping)