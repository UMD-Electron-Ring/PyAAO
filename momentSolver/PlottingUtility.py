import numpy as np
import matplotlib.pyplot as plt




class PlottingUtility:
    '''
    Helpful plots
    '''

    def __init__(self, momObj=None):
        self.momObj = momObj


    def PlotEnv(self, momObj=None):
        '''
        Plot beam envelope 
        '''

        if momObj is None:
            momObj = self.momObj
        
        xr = np.sqrt(momObj.y[0,:] + momObj.y[1,:]) * 1e3 #  <x^2> = Q+ + Q-
        yr = np.sqrt(momObj.y[0,:] - momObj.y[1,:]) * 1e3 #  <y^2> = Q+ - Q-
        plt.figure(figsize=(6,4))
        plt.plot(momObj.z, xr, label='$<x>$')
        plt.plot(momObj.z, yr, label='$<y>$')

        sfquad = xr[0] * 0.25 / np.max(momObj.kquad)
        sfsol = xr[0] * 0.25 / np.max(momObj.ksol)
        plt.plot(momObj.z, momObj.ksol * sfsol, c='m') # plot scaled magnet plots in the same plot so we know where the magnets are in the lattice
        plt.plot(momObj.z, momObj.kquad * sfquad, c='k') # plot scaled magnet plots in the same plot so we know where the magnets are in the lattice
        plt.xlabel('Z [m]')
        plt.ylabel('Beam size [mm]')
        plt.grid(True)
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
        for elem in momObj.lattice:
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

                