import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from momentSolver import MomentSolver,OptimizationUtility,MomentSolverUtility,PlottingUtility
                
# Magnet parameters & Opt parameters
from magnetParameters import lattice
from optParameters import params,parammapping

plt.close('all')

# for plotting
csfont = {"fontname":"Times New Roman"}
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["text.usetex"] = True
# colors for plotting things
colors =  plt.rcParams["axes.prop_cycle"].by_key()["color"]
# plotting stuff
title_size = 28
label_size = 24
leg_size   = 16
tick_size = 16
GR = 1.618
COLOR = "black"
mpl.rcParams["text.color"] = COLOR
mpl.rcParams["axes.labelcolor"] = COLOR
mpl.rcParams["xtick.color"] = COLOR
mpl.rcParams["ytick.color"] = COLOR

ls = ["solid","dashed","dotted","dashdot"]
mark = ['+', '.', 'o', '*']
pnames = ["Init_Sol",
          "Opt_Sol",
          "Q_Sol",
          "P_Sol",
          "E_Sol",
          "Params",
          "FoM",
          "FoM_i"
          ]

# units for easy conversion
um = 1.0e-6; mm = 1.0e-3; us = um; ms = mm

# restrictions on the optimization. This gets applied during each gradient step
def setRestrictions(momObj, paramArray):

    # right now false, so there is no restriction. If we want to add one just make this true and put in the logic
    if False:
        # We want quad 1 & 3 to have same strengths, so in our case, have the same optimization values        
        # Thus, in our optimization parameters, quad 1 strength and quad 3 strength were parameters 0 and 2
        paramArray[0] = paramArray[2]# 7/21 I want the 1st and the 3rd to be the same

    # other examples of restrictions you can do:
    # If you wanted to restrict the position of magnets
    # say you want quad 1 and quad 2 zstart positions to be atleast 0.18meters away, we can do some logic like:
    # if False:
    #     if ( paramArray[5] * mom.latticeDefault[3]['zstart'] - mom.lattice[2]['zend'] < 0.18 ):
    #         paramArray[5] = (0.18 + mom.lattice[2]['zend']) / mom.latticeDefault[3]['zstart']
    # Logic is a bit more complicated here since the optimization parameters are scalars, we need to get the actual z positions and multiply by scalars
    # Let me know if you want to do stuff like this and I can help
    # I am leaving the below as kinda more examples
    # if False:
    #     length = 0.12 # 12 cm
    #     quad3Center = momObj.lattice[2]['zstart'] * paramArray[6] + momObj.lattice[2]['length'] / 2.0
    #     paramArray[9] = (length + quad3Center) / (momObj.latticeDefault[-1]['zstart'])
    #     #paramArray[9] = 0.8
    #     #print( "POS: " + str(quad3Center) + "|" + str(paramArray[9] * momObj.lattice[3]['zstart']) )        

    return paramArray, momObj

msu = MomentSolverUtility()
ou = OptimizationUtility(setRestrictions) # restrictions function handle gets set here
pu = PlottingUtility()

# initial values
betax,betay = 0.629, 0.0629
alphax,alphay = 0.0,0.0
emitx,emity = 53.0e-6, 5.3e-6
initCond = msu.GetInitialConditions(betax,betay,alphax,alphay,emitx,emity)

# physics settings
energy = 5e3 # [eV]
current = 2.9e-3# [Amps]
pipeRadius = 0.0 # [meters] , for image charges effect on pipe walls, zero ignores the effect

# sim parameters
zF = 1.4 # final length in z where to optimize
zInterval = (0, zF)
stepSize = 0.0001

print('Setting up initial lattice')
mom = MomentSolver(lattice, energy=energy, current=current, pipeRadius=pipeRadius, zInterval=zInterval, stepSize=stepSize)
mom.initialMoments = initCond

# lattice gets initialize with the initial magnet values from the magnetParameters file.
# if we want to then update these with initial optimization parameters, in case the opt parameters aren't all 1.0 starting, we do it like this:
mom.UpdateLattice(params = params)
pu.printLattice(mom) # nice tabular printout 

pname = r"$\varepsilon=$%.2e $\beta_x=$%.2f $E_k=$%.2e $I_b=$%.2e $\alpha_x=$%.2f $z_F=$%.2f"%(emitx,betax,energy,current,alphax,zF)

# run moment equations and plot whatever we have
mom.Run()
pu.PlotEnv(mom, title='Init. Sol. ' + pname)

###############################################################

# run moment solver. Ctrl-c to interrupt
mom, an_h, gamma_h, f_h, fp_h, df_h  = ou.runOptimization(mom, params, parammapping)

# save results
# an_h is an history of the value of optimization parameters over every iteration step
# gamma_h is the history of the gamma parameters every iteration step
# f_h is the history of the FoM every iteration step
# fp_h is the history of each component of the FoM every iteration step, this is used mostly for debugging to see how much each components contributes to optimization
# df_h is gradient used in the gradient descent optimization, every iteration step
# I usually save data like this incase I need to reload old optimization values/parameters and start a new run from those values
#np.save('run1.npy',[an_h,gamma_h,f_h,fp_h,df_h],dtype=object)

# make your own plots here if you want
if False:
    plt.figure()
    plt.plot(mom.z, mom.y[0,:]) # plots Q+ vs z
    x2 = mom.y[0,:] + mom.y[1,:] # Q+ + Q- = <x^2>
    plt.plot(mom.z, x2) # plots <x^2> vs z

# print final results 
print(an_h[-1]) # print the last best set of optimization parameters
an_h[-1],mom = setRestrictions(mom, an_h[-1]) # apply the restriction function to these last set of opt parameters
mom.UpdateLattice( params = ou.getParamObj(an_h[-1],parammapping) ) # update the lattice with the last set of opt parameters
pu.printLattice(mom) # print out the lattice values

# plot stuff again, post optimization
mom.Run()
pu.PlotEnv(mom, title='Opt. Sol. ' + pname)

plt.figure() # plot Q's
plt.plot(mom.z, mom.y[0,:]/um,linestyle='-',marker=mark[0],label=r'$Q_+$')
plt.plot(mom.z, mom.y[1,:]/um,linestyle='--',marker=mark[1],label=r'$Q_-$')
plt.plot(mom.z, mom.y[2,:]/um,linestyle='-.',marker=mark[2],label=r'$Q_x$')
plt.grid(True)
plt.xlabel(r"$z$ [m]")
plt.ylabel(r"$Q$ [mm$^2$]")
plt.legend(loc="best")
plt.tight_layout()

plt.figure() # plot P's
plt.plot(mom.z, mom.y[3,:]/um,linestyle='-',marker=mark[0],label=r'$P_+$')
plt.plot(mom.z, mom.y[4,:]/um,linestyle='--',marker=mark[1],label=r'$P_-$')
plt.plot(mom.z, mom.y[5,:]/um,linestyle='-.',marker=mark[2],label=r'$P_x$')
plt.plot(mom.z, mom.y[9,:]/um,linestyle=':',marker=mark[2],label=r'$L$')
plt.grid(True)
plt.xlabel(r"$z$ [m]")
plt.ylabel(r"$P$ and $L$ [mm mrad]")
plt.legend(loc="best")
plt.tight_layout()

plt.figure()
plt.plot(mom.z, mom.y[6,:]/um,linestyle='-',marker=mark[0],label=r'$E_+$')
plt.plot(mom.z, mom.y[7,:]/um,linestyle='--',marker=mark[1],label=r'$E_-$')
plt.plot(mom.z, mom.y[8,:]/um,linestyle='-.',marker=mark[2],label=r'$E_x$')
plt.grid(True)
plt.xlabel(r"$z$ [m]")
plt.ylabel(r"$E$ [mrad$^2$]")
plt.legend(loc="best")
plt.tight_layout()

# Plot History of Parameter Variation
pu.plotParams(an_h)

# plto Figure of Merit
pu.plotfom(f_h,fp_h)

plt.show()

# save figures (all 8 of them as .pngs)
dname = "./../../../Solenoid Vary/Ax=%.2f_Bx=%.2f/"%(alphax,betax)
if (os.path.exists(dname)==0):
    os.mkdir(dname)
    
    
fig_nums = plt.get_fignums()
figs = [plt.figure(n) for n in fig_nums]
for i,fig in enumerate(figs):
    fig.savefig(dname + pnames[i] + ".png", dpi=200,transparent=True)

# save data
np.save(dname + "Parameters.npy",np.array([an_h,gamma_h,f_h,fp_h,df_h],dtype=object))

# EOF