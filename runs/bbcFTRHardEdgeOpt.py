import numpy as np
import matplotlib.pyplot as plt
from momentSolver.MomentSolver import MomentSolver,OptimizationUtility,MomentSolverUtility
from momentSolver.PlottingUtility import PlottingUtility
# Magnet parameters & Opt parameters
from systems.bbcFTRHardEdge.magnetParameters import lattice
from systems.bbcFTRHardEdge.optParameters import params,parammapping

msu = MomentSolverUtility()
ou = OptimizationUtility()
pu = PlottingUtility()

def setRestrictions(momObj, paramArray):

    if True:
        paramArray[-1] = np.clip(paramArray[-1],1,10)
    if False:
        for i,param in enumerate(paramArray[0:-1]):
            paramArray[i] = np.clip(param, 0, 10)

    if False:
        length = 0.12 # 12 cm
        quad3Center = momObj.lattice[2]['zstart'] * paramArray[6] + momObj.lattice[2]['length'] / 2.0
        paramArray[9] = (length + quad3Center) / (momObj.latticeDefault[-1]['zstart'])
        #paramArray[9] = 0.8
        #print( "POS: " + str(quad3Center) + "|" + str(paramArray[9] * momObj.lattice[3]['zstart']) )

    return paramArray

# FoM to use
useFoMForMatching = False

# initial values
betax,betay = 0.698,0.698
alphax,alphay = 0,0
emitx,emity = 53e-6,5.3e-6
initCond = msu.GetInitialConditions(betax,betay,alphax,alphay,emitx,emity)

# physics settings
energy = 5e3 # eV
current = 0.0e-3 # Amps
pipeRadius = 0.0 # meters , for image charges effect on pipe walls

# sim parameters
zInterval = (0, 1.36)
stepSize = 0.0001

print('Setting up initial lattice')
mom = MomentSolver(lattice, energy=energy, current=current, pipeRadius=pipeRadius, zInterval=zInterval, stepSize=stepSize)
mom.initialMoments = initCond
mom.UpdateLattice(params = params)
pu.printLattice(mom)

# plot initial stuff
mom.Run()
pu.PlotEnv(mom)
mom.zInterval = (0, 2.5)
mom.Run()
pu.PlotEnv(mom)
mom.zInterval = zInterval

# run moments and adjoint equations
print('Running Mom. Eqn.')
mom.Run()
mom.RunAdjoint(useMatchFoM=useFoMForMatching)


# get FoM
print('Starting Opt.')
f0,f0p,_ = mom.GetFoM_And_DFoM()
df0 = mom.GetDF()
gamma = f0 / np.sum( df0**2 )
print('Starting FoM: ' + str(f0))

# opt history
an_h = [ou.getParamArray(params)]
gamma_h = [gamma]
f_h = [f0]
fp_h = [f0p]
df_h = [df0]

# initial first step
an_h.append( an_h[0] - gamma_h[0] * df_h[0] )
mom.UpdateLattice( params = ou.getParamObj(an_h[-1],parammapping) )
an_h[-1] = setRestrictions(mom, an_h[-1])
mom.UpdateLattice( params = ou.getParamObj(an_h[-1],parammapping) )
mom.Run()
ftmp,fptmp,_ = mom.GetFoM_And_DFoM()
f_h.append(ftmp)
fp_h.append(fptmp)
print('FoM: ' + str(f_h[-1]))

# find the starting gamma value
while f_h[-1] >= f0:
    gamma_h.append( gamma_h[-1] / 2.0 )
    an_h.append( an_h[0] - gamma_h[-1] * df_h[0] )
    an_h[-1] = setRestrictions(mom, an_h[-1])    
    mom.UpdateLattice( params = ou.getParamObj(an_h[-1],parammapping) )
    mom.Run()
    ftmp,fptmp,_ = mom.GetFoM_And_DFoM()
    f_h.append(ftmp)
    fp_h.append(fptmp)
    print('FoM: ' + str(f_h[-1]))

# main loop 
try:
    while True:

        # step
        ii=1
        while f_h[-1] < f_h[-2]:
            print('Iterating ' + str(ii))

            # iterate
            an_h.append( an_h[-1] - gamma_h[-1] * df_h[-1] )       
            an_h[-1] = setRestrictions(mom, an_h[-1])        
            mom.UpdateLattice( params = ou.getParamObj(an_h[-1],parammapping) )
            mom.Run()
            ftmp,fptmp,_ = mom.GetFoM_And_DFoM()
            f_h.append(ftmp)
            fp_h.append(fptmp)
            print('FoM: ' + str(f_h[-1]))
            ii += 1

            # if we have done 20 steps in a row, start increasing the step size
            if( ii > 20 ):
                gamma_h.append( gamma_h[-1] * 2.0 )

        # can't step anymore, recompute adjoint
        print('Recomputing Adjoint Equations')

        # grab last good setting
        an_h.append( an_h[-2] )
        f_h.append( f_h[-2] )
        fp_h.append( fp_h[-2] )

        # calculate adjoint
        mom.UpdateLattice( params = ou.getParamObj(an_h[-1],parammapping) )
        mom.Run()
        mom.RunAdjoint(useMatchFoM=useFoMForMatching)

        # calculate df
        df_h.append( mom.GetDF() )

        # no improvement from last step, try to update gamma
        if( ii == 2 ):
            print('Updating Gamma')
            f0n = f_h[-1]
            ann = an_h[-1]

            iii = 1
            while f_h[-1] >= f0n:
                gamma_h.append( gamma_h[-1] / 2.0 )           
                an_h.append( ann - gamma_h[-1] * df_h[-1] )
                an_h[-1] = setRestrictions(mom, an_h[-1])            
                mom.UpdateLattice( params = ou.getParamObj(an_h[-1], parammapping) )
                mom.Run()
                ftmp,fptmp,_ = mom.GetFoM_And_DFoM()
                f_h.append(ftmp)
                fp_h.append(fptmp)
                print('FoM: ' + str(f_h[-1]))
                iii += 1

                if ( iii > 25):
                    break

        if ( f_h[-1] < 1e-14 ):
            break

        if ( len(f_h) > 50000 ):
            break

        if ( ii == 1 ):
            break
except KeyboardInterrupt:
    pass

print(an_h[-1])
pu.printLattice(mom)

# plot  stuff
pu.PlotEnv(mom)
mom.zInterval = (0, 2.5)
mom.Run()
pu.PlotEnv(mom)
mom.zInterval = zInterval

plt.show()