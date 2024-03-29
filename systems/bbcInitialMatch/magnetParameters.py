import numpy as np
from Magnets import QuadProfile, SolenoidProfile

###################################
# quad profile settings
g0 = 0.0361 # T / m A
l = 0.12 # meters
s0 = 0.06 # meters
d = 0.021 # meters
# solenoid profile settings
b = 5.0408
c = 0.5027
ds = 137.08
ls = 200 # cm
s0s = 100 #cm

# Quad settings
quad1 = {}
quad1['type'] = 'quad' # what type of magnet is this
quad1['zstart'] = .00425 # starting location in meters
quad1['length'] = l # length of the quadrupole
quad1['zend'] = quad1['zstart'] + quad1['length']
quad1['rotation'] = 0.0 # rotation angle of the quadrupole in radians
quad1['dbdx'] = QuadProfile( lambda s : -g0 * np.exp(-1 * (s-s0)**2 / d**2 ) ) # field profile function of the quadrupole

quad2 = {}
quad2['type'] = 'quad'
quad2['zstart'] = 0.10655
quad2['length'] = l
quad2['zend'] = quad2['zstart'] + quad2['length']
quad2['rotation'] = 0.0
quad2['dbdx'] = QuadProfile( lambda s : g0 * np.exp(-1 * (s-s0)**2 / d**2 ) )

quad3 = {}
quad3['type'] = 'quad'
quad3['zstart'] = 0.20895
quad3['length'] = l
quad3['zend'] = quad3['zstart'] + quad3['length']
quad3['rotation'] = 0.0
quad3['dbdx'] = QuadProfile( lambda s : -g0 * np.exp(-1 * (s-s0)**2 / d**2 ) )

quad4 = {}
quad4['type'] = 'quad'
quad4['zstart'] = 0.3100
quad4['length'] = l
quad4['zend'] = quad4['zstart'] + quad4['length']
quad4['rotation'] = 0.0
quad4['dbdx'] = QuadProfile( lambda s : g0 * np.exp(-1 * (s-s0)**2 / d**2 ) )

lattice = np.array([ quad1, quad2, quad3, quad4 ])
###################################