import numpy as np
from momentSolver.Magnets import QuadProfile, SolenoidProfile

###################################
# quad profile settings
# based on Santiago's technote 0222-2023-sbds

#initial values
# q1m = -0.037913 # T / m
# q2m = 0.047133 # T /m
# q3m = -0.028188 # T / m
# q4m = -0.010840 # T / m
# q1d = -0.017229 # T / m
# q2d = 0.020904 # T / m
# q3d = -0.017229 # T / m

# values = [q1m, q2m, q3m, q4m, q1d, q2d, q3d]
# for i,value in enumerate(values):
#    start=.98*value# array notationa not needed for modifying
#    stop=1.02*value
#    values[i]=round(np.random.uniform(start, stop), 6)
# for x in values:
#    print("Random number within 2% of " + str(x) + " = ")
#    start = .98 * x# need brackets otherwise function
#    stop = 1.02 * x
#    output = round(np.random.uniform(start, stop,10), 6)
#    print(str(output) + "\n")

# q1match = values[0] # T / m
# q2match = values[1] # T /m
# q3match = values[2] # T / m
# q4match = values[3] # T / m
# q1dbdx = values[4] # T / m
# q2dbdx = values[5] # T / m
# q3dbdx = values[6] # T / m
#is the solenoid an example of one of the magnet's strength that needs to be done?

q1match = -0.037658 # T / m
q2match = 0.04702 # T /m
q3match = -0.028066# T / m
q4match = -0.010854 # T / m
q1dbdx = -0.017231 # T / m
q2dbdx = 0.02105 # T / m
q3dbdx = -0.017419 # T / m
#Dr. Bernal Santiago's phd dealing with solenoid focusing 
#solenoid strength low energy
#how does it fit? Electron Cooling. 
soldb = 6.85 * 1e-4 # T
leff = 0.05164 #0.0589 # meters # length of the quadrupole
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
print("solstart",solstart)
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
sol['dbdx'] = SolenoidProfile( lambda s,I : soldb )
lattice.append(sol)

lattice = np.array(lattice)
###################################