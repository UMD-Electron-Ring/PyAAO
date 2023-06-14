import numpy as np
from momentSolver.Magnets import QuadProfile

###################################

q1match = -0.037658 # T / m
q2match = 0.04702 # T /m
q3match = -0.028066# T / m
q4match = -0.010854 # T / m
leff = 0.05164 #0.0589 # meters # length of the quadrupole
leffHalf = leff / 2.0 # meters
offset = leffHalf # since numbers are quad centers
q1matchstart = 0.2100 - offset
q2matchstart = 0.3300 - offset
q3matchstart = 0.4500 - offset
q4matchstart = 0.5700 - offset
lattice = []

# Setup matching section quads
quad = {}
quad['type'] = 'quad' # what type of magnet is this
quad['zstart'] = q1matchstart
quad['length'] = leff # length of the quadrupole
quad['zend'] = quad['zstart'] + quad['length']
quad['rotation'] = 0.0 # rotation angle of the quadrupole in radians
quad['dbdx'] = QuadProfile( lambda s : q1match ) # field profile function of the quadrupole , hard edge, so no dependence on position
lattice.append(quad)

quad = {}
quad['type'] = 'quad' # what type of magnet is this
quad['zstart'] = q2matchstart
quad['length'] = leff # length of the quadrupole
quad['zend'] = quad['zstart'] + quad['length']
quad['rotation'] = 0.0 # rotation angle of the quadrupole in radians
quad['dbdx'] = QuadProfile( lambda s : q2match ) # field profile function of the quadrupole , hard edge, so no dependence on position
lattice.append(quad)

quad = {}
quad['type'] = 'quad' # what type of magnet is this
quad['zstart'] = q3matchstart
quad['length'] = leff # length of the quadrupole
quad['zend'] = quad['zstart'] + quad['length']
quad['rotation'] = 0.0 # rotation angle of the quadrupole in radians
quad['dbdx'] = QuadProfile( lambda s : q3match ) # field profile function of the quadrupole , hard edge, so no dependence on position
lattice.append(quad)

quad = {}
quad['type'] = 'quad' # what type of magnet is this
quad['zstart'] = q4matchstart
quad['length'] = leff # length of the quadrupole
quad['zend'] = quad['zstart'] + quad['length']
quad['rotation'] = 0.0 # rotation angle of the quadrupole in radians
quad['dbdx'] = QuadProfile( lambda s : q4match ) # field profile function of the quadrupole , hard edge, so no dependence on position
lattice.append(quad)

lattice = np.array(lattice)
###################################