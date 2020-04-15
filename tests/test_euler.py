import unittest

from QGL import *
from QGL.Euler import *
from QGL.tools.matrix_tools import *
from QGL.Cliffords import C1
import QGL.config
try:
  from helpers import setup_test_lib
except:
  from .helpers import setup_test_lib

class EulerDecompositions(unittest.TestCase):

	N_test = 1000

	def setUp(self):
		pass
		#setup_test_lib()
		#self.q1 = QubitFactory('q1')

	def test_zyz_decomp(self):
		for j in range(self.N_test):
			Uh = haar_unitary(2)
			Ux = zyz_unitary(*zyz_angles(Uh))
			assert is_close(Uh, Ux)

	def test_xyx_decomp(self):
		for j in range(self.N_test):
			Uh = haar_unitary(2)
			Ux = xyx_unitary(*xyx_angles(Uh))
			assert is_close(Uh, Ux)

	def test_diatomic_decomp(self):
		for j in range(self.N_test):
			Uh = haar_unitary(2)
			Ux = diatomic_unitary(*diatomic_angles(Uh))
			assert is_close(Uh, Ux)

	def test_xyx_cliffords(self):
		for j in range(24):
			Uxyx = xyx_unitary(*xyx_angles(C1[j]))
			assert is_close(Uxyx, C1[j]), f"{j}"

	def test_diatomic_cliffords(self):
		for j in range(24):
			Ud = diatomic_unitary(*diatomic_angles(C1[j]))
			assert is_close(Ud, C1[j]), f"{j}"

if __name__ == "__main__":
	unittest.main()
