from adjointFTRBindings import *
import numpy as np
import copy 
import matplotlib.pyplot as plt

class MomentSolver:
    '''
    Solves the moment equations
    '''

    def __init__(self, lattice=None, energy=5e3, current=0.0, pipeRadius=0.0, zInterval=(0.0, 0.322), stepSize=0.0001, initialConditions=None, ):
        '''
        energy [eV]
        current [A]
        pipeRadius [m]
        zInterval [m,m]
        '''
        if initialConditions is None:
            self.initialMoments = self.GetInitialMoments()
        else:
            self.initialMoments = initialConditions
        self.energy = energy
        self.current = current
        self.pipeRadius = pipeRadius
        self.lattice = copy.deepcopy(lattice)
        self.latticeDefault = copy.deepcopy(lattice)

        self.zInterval = zInterval
        self.stepSize = stepSize

        self.MomentSolverUtility = MomentSolverUtility()
        self.rho, self.kPerv = self.MomentSolverUtility.GetParameters(self.energy, self.current)

        # c++ bindings
        self.bindings = AdjointFTR()
        
    def GetInitialMoments(self):
        '''
        return starting conditions for the moment equations (based off Santiago's values)
        can be updated whenever.
        These represent the initial conditions of the beam for the simulations.
        If we wanted to change the starting beam size/emittance, we do it here.
        '''
        Q_plus = 0.5*(2.2581**2*1e-6 + 0.2258**2*1e-6)
        Q_minus = 0.5*(2.2581**2*1e-6 - 0.2258**2*1e-6)
        Q_x = 0.0
        P_plus = 0
        P_minus = 0
        P_x = 0
        E_plus = (7.0855**2*1e-6 + 0.70855**2*1e-6)
        E_minus = (7.0855**2*1e-6 - 0.70855**2*1e-6)
        E_x = 0
        L = 0
        phi = 0
        return np.array([Q_plus,Q_minus,Q_x,P_plus,P_minus,P_x,E_plus,E_minus,E_x,L,phi])

    def UpdateLattice(self, lattice=None, params=None):
        '''
        Update lattice parameters
        I.e. the magnet settings 
        '''

        # If we are sending in a new lattice, just update with the existing one
        if ( type(lattice) is np.ndarray ):
            self.lattice = copy.deepcopy(lattice)
            self.latticeDefault = copy.deepcopy(lattice)

        self.optParams = params

        # if we are sending in new lattice parameters, make a new lattice
        if ( type(params) is np.ndarray ):

            allPossibleOptParams = ['zstart', 'rotation', 'dbdx']
            # loop through the lattice
            for i,elem in enumerate(self.lattice):
                optParams = self.optParams[i].keys()
                # check all possible opt parameters for the element
                for optParam in allPossibleOptParams:
                    # if an opt param is in params, update it
                    if ( optParam in optParams ):

                        # special opt param
                        if ( optParam == 'dbdx' ):
                            self.lattice[i][optParam].SetScaleFactor( params[i][optParam] )
                        else:
                            self.lattice[i][optParam] = self.latticeDefault[i][optParam] * params[i][optParam]

                        # update opt param dependent values
                        if ( optParam == 'zstart'):
                            # update non opt params
                            self.lattice[i]['zend'] = self.lattice[i]['zstart'] + self.latticeDefault[i]['length']

    def Run(self, verbose=False):
        '''
        Run the moment equations over a given interval with a given step size
        Solves the following differential equations:
        dQ/dz = P
        dP/dz = E + O * Q
        dE/dz = O * P + N*L
        dL/dz = -N * Q
        See notes 1_14_21 page 3
        '''

        # run forward integration of moment equations
        z = np.arange(self.zInterval[0], self.zInterval[1] + self.stepSize, self.stepSize) # all steps

        odefunc = lambda z,Y : self.OdeMoments(z,Y,self.lattice)
        y,ksol,kquad,Omat,Nmat = self.MomentSolverUtility.ode3(odefunc,self.zInterval[0], self.stepSize, self.zInterval[1], self.initialMoments, verbose=verbose)

        # grab solutions
        self.z = z
        self.y = y
        self.ksol = ksol
        self.kquad = kquad
        self.Omat = Omat
        self.Nmat = Nmat
        
        return z,y,ksol,kquad

    def RunAdjoint(self, verbose=False, useMatchFoM=False,usePeriodicFoM=False):
        '''
        Run the adjoint moment equations backwards
        Similar to just running the regular moment equations forward
        See notes 1_14_21 page 12, middle of page set of equations
        '''

        # run forward integration of moment equations
        z = np.arange(self.zInterval[1],self.zInterval[0] - self.stepSize,-self.stepSize) # all steps

        # setup initial conditions
        if ( useMatchFoM == True ):
            _, _, self.initialMomentsAdjoint = self.GetFoM_And_DFoM_Match()
        elif ( usePeriodicFoM == True ):
            _, _, self.initialMomentsAdjoint = self.GetFoM_And_DFoM_Periodic()
        else:
            _, _, self.initialMomentsAdjoint = self.GetFoM_And_DFoM()
        initialMoments = np.concatenate((self.initialMomentsAdjoint, np.array([self.y[-1,-1]]), self.y[:,-1]))

        odefunc = lambda z,Y : self.OdeMomentsAdjoint(z,Y,self.lattice)
        y,ksol,kquad,_,_ = self.MomentSolverUtility.ode3(odefunc,self.zInterval[1], -self.stepSize, self.zInterval[0], initialMoments, verbose=verbose)

        self.zAdj = np.flip(z)
        self.yAdj = np.flip(y, 1)

        return self.zAdj, self.yAdj, ksol, kquad    

    def OdeMoments(self, z, Y, lattice):
        '''
        Main function that solves Tom's moment equations
        
        Y(1) - Q+ [L]
        Y(2) - Q- [L]
        Y(3) - Qx [L]
        Y(4) - P+ [1]
        Y(5) - P- [1]
        Y(6) - Px [1]
        Y(7) - E+ [1/L]
        Y(8) - E- [1/L]
        Y(9) - Ex [1/L]
        Y(10) - L [L]
        See notes 1_14_21 page 3
        '''

        # grab the space charge, pipe radius, rigidity etc...
        k_perv = self.kPerv
        r_pipe = self.pipeRadius
        rho = self.rho

        # if we are in a magnetic field, grab the fields
        k_sol, k_quad, psi = self.MomentSolverUtility.GetMagnetFields(z, lattice, rho)
                
        # use bindings to calculate O and N matrices from C++ code
        O_mat,N_mat = self.bindings.getONmats(k_perv, k_sol, k_quad, psi, r_pipe, Y)
        
        # System of 10 equations
        dydt = np.array([
        # dQ/dz
        Y[3],
        Y[4],
        Y[5],
        # dP/dz
        Y[6] + np.matmul( O_mat[0,:], np.reshape(Y[0:3],(3,1)) )[0],
        Y[7] + np.matmul( O_mat[1,:], np.reshape(Y[0:3],(3,1)) )[0],
        Y[8] + np.matmul( O_mat[2,:], np.reshape(Y[0:3],(3,1)) )[0],
        # dE/dz
        np.matmul( O_mat[0,:], np.reshape(Y[3:6],(3,1)) )[0] + N_mat[0,0]*Y[9],
        np.matmul( O_mat[1,:], np.reshape(Y[3:6],(3,1)) )[0] + N_mat[1,0]*Y[9],
        np.matmul( O_mat[2,:], np.reshape(Y[3:6],(3,1)) )[0] + N_mat[2,0]*Y[9],
        # dL/dz
        -1.0*np.matmul( N_mat.T, Y[0:3] )[0],
        # dphi/dz
        -1.0*k_sol/2.0
        ])
        
        return dydt,k_sol,k_quad,O_mat,N_mat

    def OdeMomentsAdjoint(self, z, Yt, lattice):  
        '''
        Main function to solve the adjoint equations (backwards)
        
        Here we have the original 11 moments + the 11 adjoint moments for a total of a 22 variable ode solve
        See notes 1_14_21 page 12, middle of page set of equations
        '''
        Y = Yt[0:11] # adjoint variables
        Y2 = Yt[11:] # moment variables
        
        # grab the space charge, pipe radius, rigidity etc...
        k_perv = self.kPerv
        r_pipe = self.pipeRadius
        rho = self.rho
        
        # if we are in a magnetic field, grab the fields        
        k_sol, k_quad, psi = self.MomentSolverUtility.GetMagnetFields(z, lattice, rho)
                
        # use bindings to calculate O and N matrices from C++ code
        O_mat,N_mat = self.bindings.getONmats(k_perv, k_sol, k_quad, psi, r_pipe, Y2)
        
        # use bindings to calculate special matrix due to space charge variations from C++ code
        Mq,Mp,Mn = self.bindings.getSCVM(k_perv, Y2)

        # adjoint dE_dot / dz special
        dedot = np.matmul( Y[3:6], Mq ) + Y[9] * np.matmul( Y2[0:3], Mn ) - np.matmul( Y[0:3], Mp ) - np.matmul( Y[0:3], Mn ) * Y2[9]

        # System of 20 equations
        dydt = np.array([
        # adjoint dQ/dz
        Y[3],
        Y[4],
        Y[5],
        # adjoint dP/dz
        Y[6] + np.matmul( O_mat[0,:], np.reshape(Y[0:3],(3,1)) )[0],
        Y[7] + np.matmul( O_mat[1,:], np.reshape(Y[0:3],(3,1)) )[0],
        Y[8] + np.matmul( O_mat[2,:], np.reshape(Y[0:3],(3,1)) )[0],
        # adjoint dE/dz
        np.matmul( O_mat[0,:], np.reshape(Y[3:6],(3,1)) )[0] + N_mat[0,0]*Y[9] + dedot[0],
        np.matmul( O_mat[1,:], np.reshape(Y[3:6],(3,1)) )[0] + N_mat[1,0]*Y[9] + dedot[1],
        np.matmul( O_mat[2,:], np.reshape(Y[3:6],(3,1)) )[0] + N_mat[2,0]*Y[9] + dedot[2],
        # adjoint dL/dz
        -1.0*np.matmul( N_mat.T, Y[0:3] )[0],        
        # adjoint dphi/dz
        -1.0*k_sol/2.0,        
        # dQ/dz
        Y2[3],
        Y2[4],
        Y2[5],
        # dP/dz
        Y2[6] + np.matmul( O_mat[0,:], np.reshape(Y2[0:3],(3,1)) )[0],
        Y2[7] + np.matmul( O_mat[1,:], np.reshape(Y2[0:3],(3,1)) )[0],
        Y2[8] + np.matmul( O_mat[2,:], np.reshape(Y2[0:3],(3,1)) )[0],
        # dE/dz
        np.matmul( O_mat[0,:], np.reshape(Y2[3:6],(3,1)) )[0] + N_mat[0,0]*Y2[9],
        np.matmul( O_mat[1,:], np.reshape(Y2[3:6],(3,1)) )[0] + N_mat[1,0]*Y2[9],
        np.matmul( O_mat[2,:], np.reshape(Y2[3:6],(3,1)) )[0] + N_mat[2,0]*Y2[9],
        # dL/dz
        -1.0*np.matmul( N_mat.T, Y2[0:3] )[0],
        # dphi/dz
        -1.0*k_sol/2.0,
        ])

        return dydt,k_sol,k_quad,O_mat,N_mat        

    def GetCOM(self, y):
        '''
        Calculate the constant of motion
        This parameter is just a check to see if the equations are getting calculated correctly.
        If this value is not constant vs time, something is wrong.
        
        See notes 1_14_21 page 5
        '''
        L = y[9,:]
        EQ = y[6,:]*y[0,:] + y[7,:]*y[1,:] + y[8,:]*y[2,:]
        PP = y[3,:]**2 + y[4,:]**2 + y[5,:]**2
        motion = EQ + (0.5)*L**2 - (0.5)*PP
        
        return motion

    def GetFoM_And_DFoM(self, index = -1, k0 = 7):
        '''
        calculates the Figure of Merit + adjoint equation initial conditions given a set of moment values
        
        i.e.
        Given Q+,Q-,Qx,P+,P-,Px,E+,E-,Ex,L , calculate FOM:
        
        See notes 1_14_21 page 15 & 16
        '''
        
        # grab the solenoid strength
        komega = self.ksol[index]
        y = self.y # grab adjoint solutions
        e1 = 1.0 # fitting parameter 1
        e2 = 0.0 # fitting parameter 2

        f5_tmp = y[6,index] + 0.5 * komega**2 * y[0,index] - komega * y[9,index]
        f4_tmp = y[6,index] - 0.5 * komega**2 * y[0,index] + self.kPerv      

        # figure of merit broken into pieces for ease of reading
        FoM1 = 0.5 * ( y[4,index]**2 + y[5,index]**2 + y[3,index]**2 )
        FoM2 = 0.5 * (k0**2) * ( y[1,index]**2 + y[2,index]**2 )
        FoM3 = 0.5 * (k0**(-2)) * ( y[7,index]**2 + y[8,index]**2 )
        FoM4 = 0.5 * (k0**(-2)) * e1 * f4_tmp**2
        FoM5 = 0.5 * (k0**(-2)) * e2 * f5_tmp**2
            
        FoM = FoM1 + FoM2 + FoM3 + FoM4 + FoM5
        FoMp = np.array([FoM1,FoM2,FoM3,FoM4,FoM5])  

        # now calculate derivatives for adjoint variables
        # adjoint variables calculated from FoM
        dP_p = y[3,index]
        dP_m = y[4,index]
        dP_x = y[5,index]
        
        dE_p = (k0**(-2)) * e2 * (f5_tmp) * (0.5*komega**2) * (-1) + (k0**(-2)) * e1 * (f4_tmp) * (-0.5*komega**2) * (-1)
        dE_m = -k0**(2)*y[1,index]
        dE_x = -k0**(2)*y[2,index]
        
        dQ_p = (k0**(-2)) * e2 * (f5_tmp) * (-1) + (k0**(-2)) * e1 * (f4_tmp) * (-1)
        dQ_m = -k0**(-2)*y[7,index]
        dQ_x = -k0**(-2)*y[8,index]
        
        dL = (k0**(-2)) * e2 * (f5_tmp) * (-1) * komega
        
        dFoM = np.array([dQ_p,dQ_m,dQ_x,dP_p,dP_m,dP_x,dE_p,dE_m,dE_x,dL])          
            
        return FoM, FoMp, dFoM

    def GetFoM_And_DFoM_Match(self, index = -1, k0=1.0/.6489):
        '''
        FoM for a periodic match, i.e. initial conditions = final conditions
        '''
        y = self.y # grab adjoint solutions

        FoM1 = 0.5 * (k0**2) * ( (y[0,index] - y[0,0])**2 + (y[1,index] - y[1,0])**2 + (y[2,index] - y[2,0])**2 )
        FoM2 = 0.5 * ( (y[3,index] - y[3,0])**2 + (y[4,index] - y[4,0])**2 + (y[5,index] - y[5,0])**2 )
        FoM3 = 0.5 * (k0**(-2)) * ( (y[6,index] - y[6,0])**2 + (y[7,index] - y[7,0])**2 + (y[8,index] - y[8,0])**2 )
        FoM4 = 0.5 * ( (y[9,index] - y[9,0])**2)                        

        FoM = FoM1 + FoM2 + FoM3 + FoM4
        FoMp = np.array([FoM1,FoM2,FoM3,FoM4])  

        gP_p = (y[3,index] - y[3,0])
        gP_m = (y[4,index] - y[4,0])
        gP_x = (y[5,index] - y[5,0])

        gE_p = (k0**(-2)) * (y[6,index] - y[6,0])
        gE_m = (k0**(-2)) * (y[7,index] - y[7,0])
        gE_x = (k0**(-2)) * (y[8,index] - y[8,0])      

        gQ_p = (k0**2) * (y[0,index] - y[0,0])
        gQ_m = (k0**2) * (y[1,index] - y[1,0])
        gQ_x = (k0**2) * (y[2,index] - y[2,0])    

        gL = (y[9,index] - y[9,0]) 

        # now calculate derivatives for adjoint variables
        # adjoint variables calculated from FoM

        dQ_p = - gE_p
        dQ_m = - gE_m
        dQ_x = - gE_x

        dP_p = gP_p
        dP_m = gP_m
        dP_x = gP_x 

        dE_p = - gQ_p
        dE_m = - gQ_m
        dE_x = - gQ_x      

        dL = - gL         
        
        dFoM = np.array([dQ_p,dQ_m,dQ_x,dP_p,dP_m,dP_x,dE_p,dE_m,dE_x,dL])          
            
        return FoM, FoMp, dFoM

    def GetFoM_And_DFoM_Periodic(self,index=-1,k0=1.0/.6489):
        '''
        Calculates the Figure of Merit along some length in the simulation
        
        Y(1) - Q+ [L]
        Y(2) - Q- [L]
        Y(3) - Qx [L]
        Y(4) - P+ [1]
        Y(5) - P- [1]
        Y(6) - Px [1]
        Y(7) - E+ [1/L]
        Y(8) - E- [1/L]
        Y(9) - Ex [1/L]
        Y(10) - L [L]

        Returns
        -------
        returns FoM, FoMp, and dFoM

        '''
        y = self.y # grab adjoint solutions
        
        uarray = np.array([k0**2, k0**2, k0**2,
                           1.0, 1.0, 1.0,
                           1.0/k0**2, 1.0/k0**2, 1.0/k0**2, 
                           k0**2]) # array of units for multiplication of stuff
            
        FoM = np.sum((y[:-1,index] - y[:-1,0])**2*uarray)/2.0
        FoMp = np.array([FoM])  
        
        dFoM = np.array([]) # empty array for no purpose but placeholder in return method     
            
        return FoM, FoMp, dFoM

    def GetDF(self):
        '''
        Calculate Omat,Nmat perturbation matrices
        See notes 1_14_21 page 16, In order to calculate the intergral, we need to find O,N perturbation matrices.
        I.e. if I perturb a magnet, how much does O & N matrices change.
        '''
        pert = 0.01
        OmatDefault = self.Omat
        NmatDefault = self.Nmat
        paramsDefault = copy.deepcopy(self.optParams)
        params = None
        df = []
        for i,elem in enumerate(paramsDefault):
            for elemPropName in elem:
                params = copy.deepcopy(paramsDefault)
                params[i][elemPropName] = paramsDefault[i][elemPropName] + paramsDefault[i][elemPropName] * pert
                self.UpdateLattice(params = params)
                Omat,Nmat = self.CalcON(self.lattice)
                for ii in range(len(self.z)):                   
                    Omat[ii] = Omat[ii] - OmatDefault[ii]
                    Nmat[ii] = Nmat[ii] - NmatDefault[ii]  

                tmp = self.CalcInt(Omat, Nmat)
                df.append(tmp)
        return np.array(df)

    def CalcON(self, lattice):
        '''
        Calculate O and N matrices
        '''
        O_mat = np.empty(len(self.z), dtype=np.ndarray)
        N_mat = np.empty(len(self.z), dtype=np.ndarray)

        # grab the space charge, pipe radius, rigidity etc...
        k_perv = self.kPerv
        r_pipe = self.pipeRadius
        rho = self.rho
        for i,z in enumerate(self.z):
            Y = self.y[:,i]
            
            k_sol, k_quad, psi = self.MomentSolverUtility.GetMagnetFields(z, lattice, rho)
        
            # use bindings to calculate O and N matrices from C++ files
            O_mat[i],N_mat[i] = self.bindings.getONmats(k_perv, k_sol, k_quad, psi, r_pipe, Y)

        return O_mat, N_mat

    def CalcInt(self, Omat, Nmat):
        '''
        Calculate adjoint integral
        See notes 1_14_21 page 16
        '''
        int1 = np.zeros(len(self.z))
        int2 = np.zeros(len(self.z))
        int3 = np.zeros(len(self.z))
        int4 = np.zeros(len(self.z))

        y = self.y
        yadj = self.yAdj
        xx = []
        for i in range(len(self.z)):
            int1[i] = np.matmul( yadj[3:6,i], np.matmul(Omat[i], y[0:3,i]) )
            int2[i] = yadj[9,i] * np.matmul( y[0:3,i], Nmat[i] )
            int3[i] = -1 * np.matmul( yadj[0:3,i], np.matmul(Omat[i], y[3:6,i]) )
            int4[i] = -1 * np.matmul( yadj[0:3,i], Nmat[i] * y[9,i] )

        import matplotlib.pyplot as plt
        return np.trapz(int1+int2+int3+int4, self.z)

    def PrintLattice(self):
        '''
        Nicely print out lattice parameters
        '''
        pass

class MomentSolverUtility:
    '''
    Helper class with some utilities used within the MomentSolver class.
    '''

    def __init__(self):
        pass

    def GetParameters(self, energy=5e3,current=0.0):
        '''
        Grab beam ridgidity
        '''
        # Parameters
        e         = 1.60217733E-19 #C
        m         = 9.1093897E-31 #kg
        Energy    = energy # eV
        c         = 2.997924E8 # m/s

        gamma     = 1+((Energy)/(510998.9461))
        beta      = np.sqrt((gamma*gamma)-1)/gamma
        v         = beta*c
        bg        = beta*gamma
        rho       = bg*c*(m/e) 
        
        k_perv = (1.0/(4.0*np.pi))*(c*377.0*current) / (m*v**3*gamma**3/e);   
        
        return rho,k_perv

    def GetMagnetFields(self, z, lattice, rho):
        '''
        Helper function that grabs magnet field values at the given "z" position.
        '''
        # are we in magnets? get strengths
        k_sol = 0.0
        k_quad = 0.0
        psi = 0.0        
        for elem in lattice:
            # find solenoid
            if ( elem['type'] == 'solenoid' ):
                # are we within the field?
                if ( z >= elem['zstart'] and z < elem['zend'] ):
                    # add the field
                    k_sol = k_sol + (elem['dbdx'].GetValue( z - elem['zstart'] ) / rho)

            # find quadrupole
            if ( elem['type'] == 'quad' ):
                # are we within the field?
                if ( z >= elem['zstart'] and z < elem['zend'] ):
                    # add the field
                    k_quad = k_quad + (elem['dbdx'].GetValue( z - elem['zstart'] ) / rho)
                    psi = elem['rotation']

        return k_sol, k_quad, psi        

    def ode3(self, F,t0,h,tfinal,y0,verbose=False):
        '''
        third order classical Runge-Kutta ODE solver
        '''
        y = y0
        tsteps = np.arange(t0, tfinal, h)           

        # for extra params
        ksol = np.zeros(len(tsteps)+1)
        kquad = np.zeros(len(tsteps)+1)
        Omat = np.empty(len(tsteps)+1, dtype=np.ndarray)
        Nmat = np.empty(len(tsteps)+1, dtype=np.ndarray)

        _,sol,quad,Om,Nm = F(t0,y0)
        ksol[0] = sol
        kquad[0] = quad
        Omat[0] = Om
        Nmat[0] = Nm

        yout = np.zeros((len(y0),len(tsteps)+1))
        yout[:,0] = y

        if verbose:
            N = len(tsteps)
            NN = int(N/10.0)
            for i,t in enumerate(tsteps):
                t1,sol,quad,Om,Nm = F(t,y)
                s1 = h*t1
                t2,_,_,_,_ = F(t+h/2.0, y+s1/2.0)
                s2 = h*t2
                t3,_,_,_,_ = F(t+h, y-s1+2.0*s2)
                s3 = h*t3
                y = y + (s1 + 4.0*s2 + s3)/6.0
                yout[:,i+1] = y
                ksol[i+1] = sol
                kquad[i+1] = quad
                Omat[i+1] = Om
                Nmat[i+1] = Nm
                if i % NN == 0:
                    zstr = ' z = ' + str(round(t,4))[0:6]
                    print('{:<15}'.format(zstr) +  ' | ' + str(int(1.0*i/NN*10)) + '%')

            zstr = ' z = ' + str(round(t+h,4))[0:6]
            print('{:<15}'.format(zstr) +  ' | 100%')            
        else:
            for i,t in enumerate(tsteps):
                t1,sol,quad,Om,Nm = F(t,y)
                s1 = h*t1
                t2,_,_,_,_ = F(t+h/2.0, y+s1/2.0)
                s2 = h*t2
                t3,_,_,_,_ = F(t+h, y-s1+2.0*s2)
                s3 = h*t3
                y = y + (s1 + 4.0*s2 + s3)/6.0
                yout[:,i+1] = y
                ksol[i+1] = sol
                kquad[i+1] = quad
                Omat[i+1] = Om
                Nmat[i+1] = Nm             
                
        return yout,ksol,kquad,Omat,Nmat

    # initial beam conditions from courant snyder parameters
    def GetInitialConditions(self, betax, betay, alphax, alphay, emitx, emity):
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
    
    # initial beam conditions from moments (usually given from Warp)
    def GetInitialConditionsFromWarp(self, xsq, xprimesq, ysq, yprimesq, xy, xxprime, yyprime, xyprime, yxprime, xprimeyprime):
        Q_plus = 0.5*(xsq + ysq)
        Q_minus = 0.5*(xsq - ysq)
        Q_x = xy
        P_plus = xxprime + yyprime
        P_minus = xxprime - yyprime
        P_x = yxprime + xyprime
        E_plus = xprimesq + yprimesq
        E_minus = xprimesq - yprimesq
        E_x = 2.0 * xprimeyprime
        L = xyprime - yxprime
        phi = 0
        return np.array([Q_plus,Q_minus,Q_x,P_plus,P_minus,P_x,E_plus,E_minus,E_x,L,phi])        

class OptimizationUtility:
    '''
    Helper class for running optimizations
    '''

    def __init__(self, restFunc=None):
        self.restFunc = restFunc

    # get optimization parameter array from an param object
    def getParamArray(self, paramObj):
        arr = []
        for i,elem in enumerate(paramObj):
            for elemName in elem:
                arr.append(paramObj[i][elemName])
        return np.array(arr)

    # get optimization parameter object from array and param mapping
    def getParamObj(self, paramArray, paramMapping):
        params = []
        cc = 0
        for ii,elem in enumerate(paramMapping):
            tmp = {}
            for jj,param in enumerate(elem):
                tmp[param] = paramArray[cc]
                cc += 1
            params.append(tmp)
        return np.array(params)
    
    # take a step in the optimization
    def takeStep(self, mom, an0, gamma0, df0, pmapping):
        anNew = an0 - gamma0 * df0
        if( self.restFunc ):
            anNew, mom = self.restFunc(mom, anNew)
        mom.UpdateLattice( params = self.getParamObj(anNew,pmapping) )
        mom.Run()         
        ftmp,fptmp,_ = mom.GetFoM_And_DFoM_Periodic() ###################### FOM call here
        return mom, anNew, ftmp, fptmp    

    def runOptimization(self, mom: MomentSolver, params0, paramMapping, maxSteps=50000):
        # run moments and adjoint equations
        print('Running Mom. Eqn.')
        mom.Run()
        mom.RunAdjoint(usePeriodicFoM=True) ###################### FOM call here. If using new FoM, go into this function and update in there

        # get FoM
        print('Starting Opt.')
        f0,f0p,_ = mom.GetFoM_And_DFoM_Periodic() ###################### FOM call here
        df0 = mom.GetDF()
        gamma = f0 / np.sum( df0**2 )
        print('Starting FoM: ' + str(f0))

        # opt history
        an_h = [self.getParamArray(params0)]
        gamma_h = [gamma]
        f_h = [f0]
        fp_h = [f0p]
        df_h = [df0]

        # initial first step
        mom, antmp, ftmp, fptmp = self.takeStep(mom, an_h[0], gamma_h[0], df_h[0], paramMapping)
        an_h.append(antmp)
        f_h.append(ftmp)
        fp_h.append(fptmp)
        print('FoM: ' + str(f_h[-1]))
        try:
        # find the starting gamma value
            while f_h[-1] >= f0:
                gamma_h.append( gamma_h[-1] / 2.0 )
                mom, antmp, ftmp, fptmp = self.takeStep(mom, an_h[0], gamma_h[-1], df_h[0], paramMapping)
                an_h.append(antmp)
                f_h.append(ftmp)
                fp_h.append(fptmp)
                print('FoM: ' + str(f_h[-1]))    

            # main loop 
            while True:

                # step
                ii=1
                while f_h[-1] < f_h[-2]:
                    print('Iterating ' + str(ii))

                    # iterate
                    mom, antmp, ftmp, fptmp = self.takeStep(mom, an_h[-1], gamma_h[-1], df_h[-1], paramMapping)
                    an_h.append(antmp)
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
                mom.UpdateLattice( params = self.getParamObj(an_h[-1],paramMapping) )
                mom.Run()
                mom.RunAdjoint(useMatchFoM=True)

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
                        mom, antmp, ftmp, fptmp = self.takeStep(mom, ann, gamma_h[-1], df_h[-1], paramMapping)
                        an_h.append(antmp)
                        f_h.append(ftmp)
                        fp_h.append(fptmp)
                        print('FoM: ' + str(f_h[-1]))  
                        iii += 1

                        if ( iii > 25):
                            break

                if ( f_h[-1] < 1e-14 ):
                    break

                if ( len(f_h) > maxSteps ):
                    break

                if ( ii == 1 ):
                    break
        except KeyboardInterrupt:
            pass  

        # opt history
        return mom, an_h, gamma_h, f_h, fp_h, df_h  

class PlottingUtility:
    '''
    Helpful plots
    '''

    def __init__(self, momObj=None):
        self.momObj = momObj

    def PlotEnv(self, momObj=None, title=None):
        '''
        Plot beam envelope 
        '''

        if momObj is None:
            momObj = self.momObj
        
        xr = np.sqrt(momObj.y[0,:] + momObj.y[1,:]) * 1e3 #  <x^2> = Q+ + Q-
        yr = np.sqrt(momObj.y[0,:] - momObj.y[1,:]) * 1e3 #  <y^2> = Q+ - Q-
        plt.figure(figsize=(6,4))
        plt.plot(momObj.z, xr, label=r'$\sqrt{\langle x^2 \rangle}$')
        plt.plot(momObj.z, yr, '--', label=r'$\sqrt{\langle y^2 \rangle}$')

        sfquad = xr[0] * 0.25 / np.max(momObj.kquad)
        try:
            sfsol = xr[0] * 0.25 / np.max(momObj.ksol)
            plt.plot(momObj.z, momObj.ksol * sfsol, c='m') # plot scaled magnet plots in the same plot so we know where the magnets are in the lattice
        except:
            pass
        plt.plot(momObj.z, momObj.kquad * sfquad, c='k') # plot scaled magnet plots in the same plot so we know where the magnets are in the lattice
        plt.xlabel('Z [m]')
        plt.ylabel('Beam size [mm]')
        plt.grid(True)
        plt.title(title)
        plt.legend()

    def plotParams(self, params):
        '''
        Plot optimization parameters vs steps
        '''

        aa = np.zeros(( len(params),len(params[0]) ))
        for i,param in enumerate(params):
            aa[i,:] = param
        
        plt.figure()
        plt.plot(aa)

    def plotfom(self, fomh, fomph=None):
        '''
        Plot optimization parameters vs steps
        '''

        ff = np.zeros(( len(fomh),1 ))
        for i,val in enumerate(fomh):
            ff[i,:] = val
        
        plt.figure()
        plt.plot(np.log10(ff),label='FoM')

        if fomph is not None:
            ff = np.zeros(( len(fomph),len(fomph[0]) ))
            for i,val in enumerate(fomph):
                ff[i,:] = val
            plt.plot(np.log10(ff),label='FoMp')

    def printParameters(self, momObj=None, params=None):
        '''
        Print optimization parameters
        '''

        if momObj is None:
            momObj = self.momObj
        if params is None:
            params = momObj.optParams        

    def printLattice(self, momObj=None, lattice=None):
        '''
        Print lattice
        '''

        try:
            from tabulate import tabulate
        except:
            print('Please install the tabulate package for this function')
            return
        
        if momObj is None:
            momObj = self.momObj
        if lattice is None:
            lattice = momObj.lattice

        table = []
        for elem in lattice:
            entry = [
                elem['type'],
                elem['zstart'],
                (elem['zstart'] + elem['zend']) * 0.5,
                elem['zend'],
                elem['dbdx'].GetValue(0),
                elem['rotation']
            ]
            table.append(entry)

        print(tabulate(table , headers=["Element", 'zStart [m]', 'zCenter [m]', 'zEnd [m]', 'dB/dx [T/m]', 'angle [rad]']))