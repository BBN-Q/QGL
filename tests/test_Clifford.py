import unittest

from QGL import *
from QGL.tools.matrix_tools import *
from QGL.tools.clifford_tools import *
import QGL.config
try:
  from helpers import setup_test_lib
except:
  from .helpers import setup_test_lib

class EulerDecompositions(unittest.TestCase):

    def setUp(self):
        self.N1 = 24 #number of single qubit Cliffords 
        self.N2 = 11520 #number of two qubit Cliffords
        self.N_test_2 = 30 #number of two qubit Cliffords to test

    def test_n_cliffords(self):
        assert len(C1Seqs) == 24 
        assert len(C2Seqs) == 11520

    def test_multiply(self):
        for j, k in product(range(self.N1), range(self.N1)):
            m = C1[clifford_multiply(j, k)]
            mtemp = (C1[k]@C1[j]).transpose().conj()
            assert np.isclose(np.abs((mtemp@m).trace()), 2.0)

    def test_inverse(self):
        for j in range(self.N1):
            inv = C1[inverse_clifford(C1[j])]
            assert is_close(inv@C1[j], pI)
        for j in np.random.choice(range(11520), self.N_test_2):
            C = clifford_mat(j, 2)
            Ci = clifford_mat(inverse_clifford(C), 2)
            assert np.isclose(np.abs((Ci@C).trace()), 4.0)     

if __name__ == "__main__":
    unittest.main()
