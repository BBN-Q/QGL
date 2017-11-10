import unittest

from QGL import *
from QGL.PulseSequencer import *
import QGL.config
try:
  from helpers import setup_test_lib
except:
  from .helpers import setup_test_lib

class PulseTypes(unittest.TestCase):
    def setUp(self):
        setup_test_lib()
        self.q1 = QubitFactory('q1')
        self.q2 = QubitFactory('q2')
        self.q3 = QubitFactory('q3')
        self.q4 = QubitFactory('q4')

    @unittest.skip("Type promition for CNOT(q1, q2) * X(q3) gives PulseBlock not CompoundGate. Looking into this issue.")
    def test_promotion_rules(self):
        q1, q2, q3, q4 = self.q1, self.q2, self.q3, self.q4

        assert( type(X(q1)) == Pulse )
        assert( type(X(q1) + Y(q1)) == CompositePulse )
        assert( type(X(q1) * X(q2)) == PulseBlock )
        assert( type((X(q1) + Y(q1)) * X(q2)) == PulseBlock )

        assert( type(CNOT_CR(q1, q2) * X(q3)) == CompoundGate )
        assert( type(X(q3) * CNOT_CR(q1, q2)) == CompoundGate )
        assert( type(CNOT_CR(q1, q2) * CNOT_CR(q3, q4)) == CompoundGate )
