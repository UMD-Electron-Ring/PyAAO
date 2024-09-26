import numpy as np
from Magnets import QuadProfile, SolenoidProfile

###################################
# quad profile settings
# Optimized Settings using FoM GetFoM_And_DFoM_Periodic()
q1dbdx = -0.0239769              # T / m
q2dbdx = 0.033353 # T / m
q3dbdx = -0.0110646 # T / m
leff = 0.05164 # meters
leffHalf = leff / 2.0 # meters
q1start = 0.11805
q2start = 0.323234
q3start = 0.510523

lattice = []

# setup FTR section
# Quad settings
quad = {}
quad['type'] = 'quad' # what type of magnet is this
quad['zstart'] = q1start
quad['length'] = leff # length of the quadrupole
quad['zend'] = quad['zstart'] + quad['length']
# quad['rotation'] = 45*np.pi/180 # rotation angle of the quadrupole in radians
quad['rotation'] = 0.0
quad['dbdx'] = QuadProfile( lambda s : q1dbdx ) # field profile function of the quadrupole , hard edge, so no dependence on position
lattice.append(quad)

quad = {}
quad['type'] = 'quad'
quad['zstart'] = q2start
quad['length'] = leff
quad['zend'] = quad['zstart'] + quad['length']
# quad['rotation'] = 45*np.pi/180
quad['rotation'] = 0.0
quad['dbdx'] = QuadProfile( lambda s : q2dbdx )
lattice.append(quad)

quad = {}
quad['type'] = 'quad'
quad['zstart'] = q3start
quad['length'] = leff
quad['zend'] = quad['zstart'] + quad['length']
# quad['rotation'] = 45*np.pi/180
quad['rotation'] = 0.0
quad['dbdx'] = QuadProfile( lambda s : q3dbdx )
lattice.append(quad)

# quad = {}
# quad['type'] = 'quad'
# quad['zstart'] = q4start
# quad['length'] = leff
# quad['zend'] = quad['zstart'] + quad['length']
# # quad['rotation'] = 45*np.pi/180
# quad['rotation'] = 0.0
# quad['dbdx'] = QuadProfile( lambda s : q4dbdx )
# lattice.append(quad)

'''
sol = {}
sol['type'] = 'solenoid'
sol['zstart'] = solstart
sol['length'] = leffSolenoid # m
sol['zend'] = sol['zstart'] + sol['length']
sol['rotation'] = 0.0 # doesn't do anything
sol['dbdx'] = SolenoidProfile( lambda s,I : soldb )
lattice.append(sol)
'''


lattice = np.array(lattice)
###################################