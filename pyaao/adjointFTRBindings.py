#
# This is the python interface/wrapper/bindings to talk to the cpp code for certain calculation intense adjoint equations
#
#
#

import ctypes
import numpy as np
import pathlib
import os

# Get the directory where this file is located
cdir = pathlib.Path(__file__).parent.resolve()

# Load the shared library
if os.name == 'nt':  # Windows
	lib = ctypes.WinDLL(cdir.parent / 'bindings' / 'AdjointFTR.dll')
else:  # Linux
	lib = ctypes.CDLL(cdir.parent / 'bindings' / 'libAdjointFTR.so')

# Define function signatures
lib.AdjointFTR_new.argtypes = []
lib.AdjointFTR_new.restype = ctypes.c_void_p

lib.AdjointFTR_destroy.argtypes = [ctypes.c_void_p]
lib.AdjointFTR_destroy.restype = None

lib.AdjointFTR_getSCVM.argtypes = [
	ctypes.c_void_p,
	np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
	ctypes.c_double,
	ctypes.c_double,
	ctypes.c_double,
	ctypes.c_double,
	ctypes.c_double,
	ctypes.c_double,
	ctypes.c_double
]
lib.AdjointFTR_getSCVM.restype = None

lib.AdjointFTR_getONmats.argtypes = [
	ctypes.c_void_p,
	np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
	ctypes.c_double,
	ctypes.c_double,
	ctypes.c_double,
	ctypes.c_double,
	ctypes.c_double,
	ctypes.c_double,
	ctypes.c_double,
	ctypes.c_double,
	ctypes.c_double
]
lib.AdjointFTR_getONmats.restype = None


class AdjointFTR:
	'''
	Python binding class for C++ AdjointFTR implementation.

	This class provides Python access to computationally intensive adjoint
	equation calculations implemented in C++.

	The object can be used as a context manager:
		with AdjointFTR() as ftr:
			Mq, Mp, Mn = ftr.getSCVM(kPerv, Y)
	'''

	def __init__(self):
		'''Initialize the C++ AdjointFTR object.'''
		self.obj = lib.AdjointFTR_new()
		if not self.obj:
			raise RuntimeError("Failed to create AdjointFTR object")

	def __del__(self):
		'''Destructor - ensures C++ object is properly cleaned up.'''
		if hasattr(self, 'obj') and self.obj:
			lib.AdjointFTR_destroy(self.obj)
			self.obj = None

	def __enter__(self):
		'''Context manager entry.'''
		return self

	def __exit__(self, exc_type, exc_val, exc_tb):
		'''Context manager exit - cleanup resources.'''
		self.__del__()
		return False

	def getSCVM(self, kPerv, Y):
		'''
		Calculate space charge variation matrices given perveance and moments.

		Parameters
		----------
		kPerv : float
			Perveance parameter
		Y : array-like
			Moment vector (must have at least 6 elements)

		Returns
		-------
		Mq : ndarray, shape (3, 3)
			First variation matrix
		Mp : ndarray, shape (3, 3)
			Second variation matrix
		Mn : ndarray, shape (3, 3)
			Third variation matrix
		'''
		if not self.obj:
			raise RuntimeError("AdjointFTR object has been destroyed")

		# Pre-allocate output buffer (27 doubles)
		out = np.zeros(27, dtype=np.float64)

		# Call C++ function - it will fill the pre-allocated buffer
		lib.AdjointFTR_getSCVM(
			self.obj, out, kPerv,
			float(Y[0]), float(Y[1]), float(Y[2]),
			float(Y[3]), float(Y[4]), float(Y[5])
		)

		# Reshape into the three matrices
		Mq = out[0:9].reshape((3, 3))
		Mp = out[9:18].reshape((3, 3))
		Mn = out[18:27].reshape((3, 3))

		return Mq, Mp, Mn

	def getONmats(self, kPerv, kSol, kQuad, kQuadRot, pipeRadius, Y):
		'''
		Calculate O and N matrices.

		Parameters
		----------
		kPerv : float
			Perveance parameter
		kSol : float
			Solenoid parameter
		kQuad : float
			Quadrupole parameter
		kQuadRot : float
			Quadrupole rotation parameter
		pipeRadius : float
			Pipe radius
		Y : array-like
			Moment vector (must have at least 11 elements)

		Returns
		-------
		Omat : ndarray, shape (3, 3)
			O matrix
		Nmat : ndarray, shape (3, 1)
			N matrix
		'''
		if not self.obj:
			raise RuntimeError("AdjointFTR object has been destroyed")

		# Pre-allocate output buffer (12 doubles)
		out = np.zeros(12, dtype=np.float64)

		# Call C++ function - it will fill the pre-allocated buffer
		lib.AdjointFTR_getONmats(
			self.obj, out, kPerv, kSol, kQuad, kQuadRot, pipeRadius,
			float(Y[0]), float(Y[1]), float(Y[2]), float(Y[10])
		)

		# Reshape into matrices
		Omat = out[0:9].reshape((3, 3))
		Nmat = out[9:12].reshape((3, 1))

		return Omat, Nmat