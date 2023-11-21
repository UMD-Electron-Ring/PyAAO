# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 15:12:35 2023

@author: lpocher
"""

import numpy as np
from Magnets import QuadProfile, SolenoidProfile

###################################
# quad profile settings
# based on Santiago's technote 0222-2023-sbds
# solutoin for 2.6 mA of space charge
q1dbdx = -0.025 # T / m
q2dbdx = 0.031156 # T / m
q3dbdx = -0.02 # T / m
# q2dbdx = 0.00 # T / m
# q3dbdx = 0.00 # T / m
# q1dbdx = -0.017 # T / m
# q2dbdx = 0.035 # T / m
# q3dbdx = -0.017 # T / m
# q4dbdx = 0.0 # T / m
q4dbdx = 0.025 # T / m
# soldb = 8.0*1.0e-4 # solenoid field in T for 2.6 mA solution
soldb = 7.0*1.0e-4 # solenoid field in T for 0.0 mA solution ~68 G
leff = 0.05164 # meters
leffHalf = leff / 2.0 # meters
lqact = 1.16*2.0*0.0254 # actual physical length of quadrupole
mSep = 0.12
quadmount = 21.25*0.0254 
q1mDist = 0.155 # minimum distance quad 1 center can get to slit  0.121539 / 0.172339
leffSolenoid = 1.37*1.0 # length of solenoid , doesn't have to be super long here, we are optimizing up to the entrance of a solenoid, so can be short
# q1start = 0.129
# q2start = 0.299466
# q3start = 0.493
# q4start = 0.58
# solstart = q1mDist + quadmount - lqact/2.0

q1start = 0.18
q2start = 0.32
q3start = 0.47
q4start = 0.59
solstart = q1mDist + quadmount - lqact/2.0
# solstart =  0.7 # solution with 2.6 mA of Space Charge
# solstart  = 0.58 # solution with 0.0 mA of Space Charge

lattice = []

# setup FTR section
# Quad settings
quad = {}
quad['type'] = 'quad' # what type of magnet is this
quad['zstart'] = q1start
quad['length'] = leff # length of the quadrupole
quad['zend'] = quad['zstart'] + quad['length']
quad['rotation'] = 45*np.pi/180 # rotation angle of the quadrupole in radians
# quad['rotation'] = 0.0
quad['dbdx'] = QuadProfile( lambda s : q1dbdx ) # field profile function of the quadrupole , hard edge, so no dependence on position
lattice.append(quad)

quad = {}
quad['type'] = 'quad'
quad['zstart'] = q2start
quad['length'] = leff
quad['zend'] = quad['zstart'] + quad['length']
quad['rotation'] = 45*np.pi/180
# quad['rotation'] = 0.0
quad['dbdx'] = QuadProfile( lambda s : q2dbdx )
lattice.append(quad)

quad = {}
quad['type'] = 'quad'
quad['zstart'] = q3start
quad['length'] = leff
quad['zend'] = quad['zstart'] + quad['length']
quad['rotation'] = 45*np.pi/180
# quad['rotation'] = 0.0
quad['dbdx'] = QuadProfile( lambda s : q3dbdx )
lattice.append(quad)

quad = {}
quad['type'] = 'quad'
quad['zstart'] = q4start
quad['length'] = leff
quad['zend'] = quad['zstart'] + quad['length']
quad['rotation'] = 45*np.pi/180
# quad['rotation'] = 0.0
quad['dbdx'] = QuadProfile( lambda s : q4dbdx )
lattice.append(quad)


sol = {}
sol['type'] = 'solenoid'
sol['zstart'] = solstart
sol['length'] = leffSolenoid # m
sol['zend'] = sol['zstart'] + sol['length']
sol['rotation'] = 0.0 # doesn't do anything
sol['dbdx'] = SolenoidProfile( lambda s : soldb )
lattice.append(sol)


lattice = np.array(lattice)
###################################