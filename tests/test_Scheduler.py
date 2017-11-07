import unittest

from QGL import *
try:
  from helpers import setup_test_lib
except:
  from .helpers import setup_test_lib

class SchedulerTest(unittest.TestCase):
    def setUp(self):
        setup_test_lib()
        self.q1 = QubitFactory('q1')
        self.q2 = QubitFactory('q2')
        self.q3 = QubitFactory('q3')
        self.q4 = QubitFactory('q4')

    def test_1q_ops(self):
        q1, q2, q3 = self.q1, self.q2, self.q3

        # uniform fill on first time step
        seq = [X(q1), Y(q2), Z(q3),
               X(q1), Y(q2)]
        result = schedule(seq)

        assert(result == [X(q1)*Y(q2)*Z(q3),
                          X(q1)*Y(q2)] )

        # add an extra pulse on q1
        seq = [X(q1), X(q1), Y(q2), Z(q3), X(q1), Y(q2)]
        result = schedule(seq)

        assert(result == [X(q1)*Y(q2)*Z(q3),
                          X(q1)*Y(q2),
                          X(q1)] )

        # same sequence but with a Barrier
        seq = [X(q1),
               Barrier(q1, q2, q3),
               X(q1), Y(q2), Z(q3),
               X(q1), Y(q2)]
        result = schedule(seq)

        assert(result == [X(q1),
                          X(q1)*Y(q2)*Z(q3),
                          X(q1)*Y(q2)] )

    def test_1q_composite(self):
        q1, q2, q3 = self.q1, self.q2, self.q3

        seq = [X(q1)+Y(q1), X(q2),
               Y(q1), Y(q2)]
        result = schedule(seq)

        assert(result == [(X(q1)+Y(q1))*X(q2),
                          Y(q1)*Y(q2)] )

    def test_2q_ops(self):
        q1, q2, q3, q4 = self.q1, self.q2, self.q3, self.q4

        seq = [X(q1), CNOT_simple(q2, q3),
               X(q1), Y(q2)]
        result = schedule(seq)

        assert(result == [X(q1)*CNOT_simple(q2, q3),
                          X(q1)*Y(q2)] )

        seq = [X(q1), CNOT_CR(q2, q3),
               X(q1), Y(q2)]
        result = schedule(seq)

        assert(result == [X(q1)*CNOT_CR(q2, q3),
                          X(q1)*Y(q2)] )

        seq = [CNOT_simple(q1, q2), CNOT_simple(q3, q4), X(q1), X(q2)]
        result = schedule(seq)

        assert(result == [CNOT_simple(q1, q2) * CNOT_simple(q3, q4),
                          X(q1) * X(q2)] )

    def test_controlflow(self):
        # TODO to do this properly we need ControlInstructions to be valid operands
        # in tensor products, i.e. we want qif(..) * qif(...) to be valid
        q1, q2, q3, q4 = self.q1, self.q2, self.q3, self.q4

        cond_seq = qif(1, [Y(q1),Y(q2),X90(q1)], [Z(q1),Z(q2),X90(q1)])
        seq = [X(q1),X(q2),
               CNOT(q1,q2)] + \
              cond_seq + \
              [Y90(q1),Y90(q2)]

        result = schedule(seq)

        # construct the scheduled conditional from the original to preserve labels
        cond_seq_scheduled = cond_seq[:2] + \
                             [Z(q1) * Z(q2), X90(q1)] + \
                             cond_seq[5:7] + \
                             [Y(q1) * Y(q2), X90(q1)] + \
                             cond_seq[-1:]
        assert(result == [X(q1) * X(q2),
                          CNOT(q1, q2)] + \
                         cond_seq_scheduled + \
                         [Y90(q1) * Y90(q2)])

        # test that "global" control flow injects barriers
        seq = [X(q1)] + \
              cond_seq + \
              [Y90(q1),Y90(q2)]

        result = schedule(seq)
        assert(result == [X(q1)] + \
                         cond_seq_scheduled + \
                         [Y90(q1) * Y90(q2)])

    def test_measurements(self):
        # test that measurements on a qubit mark that qubit as busy in a given
        # time slot
        q1 = self.q1

        seq = [X(q1), MEAS(q1), Y(q1), MEAS(q1)]
        result = schedule(seq)

        # should be unchanged
        assert seq == result
