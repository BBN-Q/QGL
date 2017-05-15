import unittest

from QGL import *
from QGL.Scheduler import schedule

class SchedulerTest(unittest.TestCase):
    def setUp(self):
        self.q1 = Qubit(label='q1')
        self.q2 = Qubit(label='q2')
        self.q3 = Qubit(label='q3')
        self.q4 = Qubit(label='q4')
        self.q1q2 = Edge(label='q1q2', source=self.q1, target=self.q2)
        self.q2q3 = Edge(label='q2q3', source=self.q2, target=self.q3)
        self.q3q4 = Edge(label='q3q4', source=self.q3, target=self.q4)

        ChannelLibrary.channelLib.channelDict = {
            'q1': self.q1,
            'q2': self.q2,
            'q3': self.q3,
            'q1q2': self.q1q2,
            'q2q3': self.q2q3,
            'q3q4': self.q3q4
        }
        ChannelLibrary.channelLib.build_connectivity_graph()

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

        # unroll the CNOT_CR sequence
        cnotseq = CNOT_CR(q2, q3)

        # FIXME currently fails because it schedules the 2nd X(q1) "inside" the
        # CNOT(q2, q3)
        # assert(result == [X(q1)*cnotseq[0],
        #                   cnotseq[1],
        #                   cnotseq[2],
        #                   cnotseq[3],
        #                   cnotseq[4],
        #                   X(q1)*Y(q2)] )

        seq = [CNOT_simple(q1, q2), CNOT_simple(q3, q4), X(q1), X(q2)]
        result = schedule(seq)

        assert(result == [CNOT_simple(q1, q2) * CNOT_simple(q3, q4),
                          X(q1) * X(q2)] )

    def test_controlflow(self):
        # TODO to do this properly we need ControlInstructions to be valid operands
        # in tensor products, i.e. we want qif(..) * qif(...) to be valid
        pass
