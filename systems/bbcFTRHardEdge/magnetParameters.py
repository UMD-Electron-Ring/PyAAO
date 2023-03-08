import sys, os
sys.path.append('..' + os.sep + '..') # add 2 dir up to path

import numpy as np
from momentSolver.Magnets import QuadProfile, SolenoidProfile

###################################
# quad profile settings
# based on Santiago's technote 0222-2023-sbds
q1dbdx = -0.017229 # T / m
q2dbdx = 0.020904 # T / m
q3dbdx = -0.017229 # T / m
soldb = 6.85 * 1e-4 # T
leff = 0.05164 # meters
leffHalf = leff / 2.0 # meters
leffSolenoid = 1.37*1.0 # length of solenoid , doesn't have to be super long here, we are optimizing up to the entrance of a solenoid, so can be short
offset = (0.6941 - leffHalf) - 0.00425 # meters
q1start = 0.6941 - offset

# Quad settings
quad1 = {}
quad1['type'] = 'quad' # what type of magnet is this
quad1['zstart'] = (0.6941 - leffHalf)  - offset # starting location in meters
quad1['length'] = leff # length of the quadrupole
quad1['zend'] = quad1['zstart'] + quad1['length']
quad1['rotation'] = 45*np.pi/180 # rotation angle of the quadrupole in radians
quad1['dbdx'] = QuadProfile( lambda s : q1dbdx ) # field profile function of the quadrupole , hard edge, so no dependence on position

quad2 = {}
quad2['type'] = 'quad'
quad2['zstart'] = (0.9128 - leffHalf) - offset
quad2['length'] = leff
quad2['zend'] = quad2['zstart'] + quad2['length']
quad2['rotation'] = 45*np.pi/180
quad2['dbdx'] = QuadProfile( lambda s : q2dbdx )

quad3 = {}
quad3['type'] = 'quad'
quad3['zstart'] = (1.1315 - leffHalf) - offset
quad3['length'] = leff
quad3['zend'] = quad3['zstart'] + quad3['length']
quad3['rotation'] = 45*np.pi/180
quad3['dbdx'] = QuadProfile( lambda s : q3dbdx )

sol = {}
sol['type'] = 'solenoid'
sol['zstart'] = (1.8573 - leffSolenoid/2.0) - offset
sol['length'] = leffSolenoid # m
sol['zend'] = sol['zstart'] + sol['length']
sol['rotation'] = 0.0 # doesn't do anything
sol['dbdx'] = SolenoidProfile( lambda s,I : soldb )

lattice = np.array([ quad1, quad2, quad3, sol ])
###################################