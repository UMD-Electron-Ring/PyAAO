import sys, os
sys.path.append('..' + os.sep + '..') # add 2 dir up to path

import numpy as np
from momentSolver.Magnets import QuadProfile, SolenoidProfile

###################################
# quad profile settings
# based on Santiago's technote 0222-2023-sbds
q1match = -0.037913 # T / m
q2match = 0.047133 # T /m
q3match = -0.028188 # T / m
q4match = -0.010840 # T / m
q1dbdx = -0.017229 # T / m
q2dbdx = 0.020904 # T / m
q3dbdx = -0.017229 # T / m
soldb = 6.85 * 1e-4 # T
leff = 0.05164 # meters
leffHalf = leff / 2.0 # meters
leffSolenoid = 1.37*1.0 # length of solenoid , doesn't have to be super long here, we are optimizing up to the entrance of a solenoid, so can be short
offset = (0.6941 - leffHalf) - 0.00425 # meters
q1start = 0.6941 - offset

lattice = []

# Setup matching section quads
quad = {}
quad['type'] = 'quad' # what type of magnet is this
quad['zstart'] = 0.2100
quad['length'] = leff # length of the quadrupole
quad['zend'] = quad['zstart'] + quad['length']
quad['rotation'] = 0.0 # rotation angle of the quadrupole in radians
quad['dbdx'] = QuadProfile( lambda s : q1match ) # field profile function of the quadrupole , hard edge, so no dependence on position
lattice.append(quad)

quad = {}
quad['type'] = 'quad' # what type of magnet is this
quad['zstart'] = 0.3300
quad['length'] = leff # length of the quadrupole
quad['zend'] = quad['zstart'] + quad['length']
quad['rotation'] = 0.0 # rotation angle of the quadrupole in radians
quad['dbdx'] = QuadProfile( lambda s : q2match ) # field profile function of the quadrupole , hard edge, so no dependence on position
lattice.append(quad)

quad = {}
quad['type'] = 'quad' # what type of magnet is this
quad['zstart'] = 0.4500
quad['length'] = leff # length of the quadrupole
quad['zend'] = quad['zstart'] + quad['length']
quad['rotation'] = 0.0 # rotation angle of the quadrupole in radians
quad['dbdx'] = QuadProfile( lambda s : q3match ) # field profile function of the quadrupole , hard edge, so no dependence on position
lattice.append(quad)

quad = {}
quad['type'] = 'quad' # what type of magnet is this
quad['zstart'] = 0.5700
quad['length'] = leff # length of the quadrupole
quad['zend'] = quad['zstart'] + quad['length']
quad['rotation'] = 0.0 # rotation angle of the quadrupole in radians
quad['dbdx'] = QuadProfile( lambda s : q4match ) # field profile function of the quadrupole , hard edge, so no dependence on position
lattice.append(quad)

# setup FTR section
# Quad settings
quad = {}
quad['type'] = 'quad' # what type of magnet is this
quad['zstart'] = 0.6941
quad['length'] = leff # length of the quadrupole
quad['zend'] = quad['zstart'] + quad['length']
quad['rotation'] = 45*np.pi/180 # rotation angle of the quadrupole in radians
quad['dbdx'] = QuadProfile( lambda s : q1dbdx ) # field profile function of the quadrupole , hard edge, so no dependence on position
lattice.append(quad)

quad = {}
quad['type'] = 'quad'
quad['zstart'] = 0.9128
quad['length'] = leff
quad['zend'] = quad['zstart'] + quad['length']
quad['rotation'] = 45*np.pi/180
quad['dbdx'] = QuadProfile( lambda s : q2dbdx )
lattice.append(quad)

quad = {}
quad['type'] = 'quad'
quad['zstart'] = 1.1315
quad['length'] = leff
quad['zend'] = quad['zstart'] + quad['length']
quad['rotation'] = 45*np.pi/180
quad['dbdx'] = QuadProfile( lambda s : q3dbdx )
lattice.append(quad)

sol = {}
sol['type'] = 'solenoid'
sol['zstart'] = lattice[-1]['zend']
sol['length'] = leffSolenoid # m
sol['zend'] = sol['zstart'] + sol['length']
sol['rotation'] = 0.0 # doesn't do anything
sol['dbdx'] = SolenoidProfile( lambda s,I : soldb )
lattice.append(sol)

lattice = np.array(lattice)
###################################