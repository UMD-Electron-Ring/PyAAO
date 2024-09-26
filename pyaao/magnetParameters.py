import numpy as np
from Magnets import QuadProfile, SolenoidProfile

###################################
# quad profile settings
g0 = 0.0361 # T/meter*A
l = 0.2 # meters
s0 = 0.06 # meters l/2.0
d = 0.021 # meters
leff = 0.05164 # meters
leffHalf = leff / 2.0 # meters
lqact = 1.16*2.0*0.0254 # actual physical length of quadrupole
mSep = 0.12
quadmount = 21.25*0.0254 
q1start = 0.23 - s0/2.0
q2start = 0.35 - s0/2.0
q3start = 0.51 - s0/2.0

# solenoid profile settings
b = 5.0408
c = 0.5027
ds = 137.08
ls = 200 # cm
s0s = 100 #cm
solstart = q1start + quadmount # physical start of the long solenoid
sollen = 1.37 # solenoid physical length in meter

# Quad settings
quad1 = {}
quad1['type'] = 'quad' # what type of magnet is this
quad1['zstart'] =q1start # starting location in meters
quad1['length'] = l # length of the quadrupole
quad1['zend'] = quad1['zstart'] + quad1['length']
quad1['rotation'] = 45*np.pi/180 # rotation angle of the quadrupole in radians
quad1['dbdx'] = QuadProfile( lambda s : -g0 * np.exp(-1 * (s-s0)**2 / d**2 ) ) # field profile function of the quadrupole

quad2 = {}
quad2['type'] = 'quad'
quad2['zstart'] = q2start
quad2['length'] = l
quad2['zend'] = quad2['zstart'] + quad2['length']
quad2['rotation'] = 45*np.pi/180
quad2['dbdx'] = QuadProfile( lambda s : g0 * np.exp(-1 * (s-s0)**2 / d**2 ) )

quad3 = {}
quad3['type'] = 'quad'
quad3['zstart'] = q3start
quad3['length'] = l
quad3['zend'] = quad3['zstart'] + quad3['length']
quad3['rotation'] = 45*np.pi/180
quad3['dbdx'] = QuadProfile( lambda s : -g0 * np.exp(-1 * (s-s0)**2 / d**2 ) )

# quad4 = {}
# quad4['type'] = 'quad'
# quad4['zstart'] = 0.3100
# quad4['length'] = l
# quad4['zend'] = quad4['zstart'] + quad4['length']
# quad4['rotation'] = 0.0
# quad4['dbdx'] = QuadProfile( lambda s : g0 * np.exp(-1 * (s-s0)**2 / d**2 ) )

# sol = {}
# sol['type'] = 'solenoid'
# sol['zstart'] = solstart + sollen/2.0 - ls*1.0e-2/2.0
# sol['length'] = ls * 1e-2 # m
# sol['zend'] = sol['zstart'] + sol['length']
# sol['rotation'] = 0.0 # doesn't do anything
# sol['dbdx'] = SolenoidProfile( lambda s,I : 1e-4*(7.4153*I - 0.12786)*c*( ((s*1e2-s0s)+ds/2)/(((s*1e2-s0s)+ds/2)**2+b**2)**0.5 - ((s*1e2-s0s)-ds/2)/(((s*1e2-s0s)-ds/2)**2+b**2)**0.5 ) )

# build your lattice, in this case only use 3 quads
# lattice = np.array([ quad1, quad2, quad3, sol])
lattice = np.array([ quad1, quad2, quad3])
###################################