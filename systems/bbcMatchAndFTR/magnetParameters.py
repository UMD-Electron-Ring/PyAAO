import numpy as np
from Magnets import QuadProfile, SolenoidProfile

###################################
# quad profile settings
# # based on Santiago's technote 0222-2023-sbds
# q1match = -0.037913 # T / m
# q2match = 0.047133 # T /m
# q3match = -0.028188 # T / m
# q4match = -0.010840 # T / m
# # q1match,q2match,q3match,q4match = 0,0,0,0
# q1dbdx = -0.017229 # T / m
# q2dbdx = 0.020904 # T / m
# q3dbdx = -0.017229 # T / m
# soldb = 6.85 * 1e-4 # T
# leff = 0.05164 # meters
# leffHalf = leff / 2.0 # meters
# leffSolenoid = 1.37*1.0 # length of solenoid , doesn't have to be super long here, we are optimizing up to the entrance of a solenoid, so can be short
# offset = leffHalf # since numbers are quad centers
# q1matchstart = 0.2100 - offset
# q2matchstart = 0.3300 - offset
# q3matchstart = 0.4500 - offset
# q4matchstart = 0.5700 - offset
# q1start = 0.6941 - offset
# q2start = 0.9128 - offset
# q3start = 1.1315 - offset
# solstart = q3start + leff + 0.18 # 18 centermeter offset

# Aita's optimized values
q1match = -0.0440226 # T / m
q2match = 0.0515796 # T /m
q3match = -0.0269442 # T / m
q4match = -0.0112155 # T / m
# q1match,q2match,q3match,q4match = 0,0,0,0
q1dbdx = -0.0174275 # T / m
q2dbdx = 0.0250344 # T / m
q3dbdx = -0.0174275 # T / m
soldb = 6.85 * 1e-4 # T
leff = 0.05164 # meters
leffHalf = leff / 2.0 # meters
leffSolenoid = 1.37*1.0 # length of solenoid , doesn't have to be super long here, we are optimizing up to the entrance of a solenoid, so can be short
offset = leffHalf # since numbers are quad centers
q1matchstart = 0.2100 - offset
q2matchstart = 0.3300 - offset
q3matchstart = 0.4500 - offset
q4matchstart = 0.5700 - offset
q1start = 0.6941 - offset
q2start = 0.9128 - offset
q3start = 1.1315 - offset
solstart = q3start + leff + 0.18 # 18 centermeter offset



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

# setup FTR section
# Quad settings
quad = {}
quad['type'] = 'quad' # what type of magnet is this
quad['zstart'] = q1start
quad['length'] = leff # length of the quadrupole
quad['zend'] = quad['zstart'] + quad['length']
quad['rotation'] = 45*np.pi/180 # rotation angle of the quadrupole in radians
quad['dbdx'] = QuadProfile( lambda s : q1dbdx ) # field profile function of the quadrupole , hard edge, so no dependence on position
lattice.append(quad)

quad = {}
quad['type'] = 'quad'
quad['zstart'] = q2start
quad['length'] = leff
quad['zend'] = quad['zstart'] + quad['length']
quad['rotation'] = 45*np.pi/180
quad['dbdx'] = QuadProfile( lambda s : q2dbdx )
lattice.append(quad)

quad = {}
quad['type'] = 'quad'
quad['zstart'] = q3start
quad['length'] = leff
quad['zend'] = quad['zstart'] + quad['length']
quad['rotation'] = 45*np.pi/180
quad['dbdx'] = QuadProfile( lambda s : q3dbdx )
lattice.append(quad)

sol = {}
sol['type'] = 'solenoid'
sol['zstart'] = solstart
sol['length'] = leffSolenoid # m
sol['zend'] = sol['zstart'] + sol['length']
sol['rotation'] = 0.0 # doesn't do anything
# sol['dbdx'] = SolenoidProfile( lambda s,I : soldb )
sol['dbdx'] = SolenoidProfile( lambda s: soldb )
lattice.append(sol)

lattice = np.array(lattice)
###################################