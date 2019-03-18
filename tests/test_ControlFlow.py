import unittest

from QGL import *
from QGL.ControlFlow import CmpEq, CmpNeq, Goto, Call, Return, LoadRepeat, Repeat
from QGL.BlockLabel import label, endlabel


class ControlFlowTest(unittest.TestCase):
    def setUp(self):
        cl = ChannelLibrary(db_resource_name=":memory:")
        self.q1 = cl.new_qubit(label='q1')
        self.q2 = cl.new_qubit(label='q2')
        self.q3 = cl.new_qubit(label='q3')
        self.q4 = cl.new_qubit(label='q4')
        cl.update_channelDict()

    def test_qif(self):
        q1 = self.q1
        seq1 = [X90(q1), Y90(q1)]
        seq2 = [X(q1), Y(q1), Z(q1)]
        label(seq1)
        label(seq2)
        assert (qif(0, seq1, seq2) == [CmpEq("m", 0), Goto(label(seq1))] + seq2 +
                [Goto(endlabel(seq1))] + seq1)

    @unittest.expectedFailure
    def test_qif_single_element(self):
        q1 = self.q1
        # just if branch
        seq = qif(0, X(q1))
        # if and else branches
        seq = qif(0, X(q1), Y(q1))

    def test_inline_qif(self):
        q1 = self.q1
        seq = [X90(q1), Y(q1), qwait(kind="CMP"), qif(0, [Id(q1)], [X(q1)]), Y(q1)]
        Compiler.compile_sequence(seq)

        seq = [X90(q1), Y(q1), qwait(kind="CMP"), qif(0, [Id(q1)], [X(q1)]), Y(q1)]
        seqs = [[Id(q1), Y(q1)], [X(q1), Y(q1)], seq]
        Compiler.compile_sequences(seqs)

    def test_qwhile(self):
        q1 = self.q1
        seq1 = [X90(q1), Y90(q1)]
        seq2 = qwhile(0, seq1)
        seq3 = [label(seq2), CmpNeq("m", 0), Goto(endlabel(seq2))] + seq1 + \
            [Goto(label(seq2)), endlabel(seq2)]
        assert (seq2 == seq3)

    def test_qdowhile(self):
        q1 = self.q1
        seq1 = [X90(q1), Y90(q1)]
        label(seq1)
        # print qdowhile(0, seq1)
        # print seq1 + [CmpEq(0), Goto(label(seq1))]
        assert (qdowhile(0, seq1) == seq1 + [CmpEq("m", 0), Goto(label(seq1))])

    def test_qcall(self):
        q1 = self.q1
        q2 = self.q2

        @qfunction
        def dummy(q):
            return [X(q), Y(q)]

        # multiple calls should return the same thing
        assert dummy(q1) == dummy(q1)
        assert dummy(q2) == dummy(q2)
        assert dummy(q1) != dummy(q2)

        # specialization lookup with label at beginning and RETURN at end
        assert ControlFlow.qfunction_specialization(dummy(
            q1).target) == [dummy(q1).target, X(q1), Y(q1), Return()]

    def test_repeat(self):
        q1 = self.q1
        seq1 = [X90(q1), Y90(q1)]
        label(seq1)
        assert (
            repeat(5, seq1) == [LoadRepeat(5)] + seq1 + [Repeat(label(seq1))])

    def test_qwait(self):
        q1 = self.q1
        seq1 = [qwait(), qwait(kind="CMP"), qwait(kind = "RAM", addr = 0)]
        assert (isinstance(seq1[0], ControlFlow.Wait))
        assert (isinstance(seq1[1], ControlFlow.LoadCmp))
        assert (isinstance(seq1[2][0], TdmInstructions.WriteAddrInstruction))
        assert (isinstance(seq1[2][1], TdmInstructions.LoadCmpVramInstruction))

    def test_qwait_err(self):
        q1 = self.q1
        with self.assertRaises(ValueError) as exc:
            seq1 = [qwait(kind='FOO')]
        exc_str = str(exc.exception)
        assert (exc_str == 'Unknown kind parameter [FOO]')

        with self.assertRaises(ValueError) as exc:
            seq1 = [qwait(kind='RAM')]
        exc_str = str(exc.exception)
        assert (exc_str == 'Please specify addr')

        # test the legal values
        seq1 = [qwait(kind='TRIG')]
        seq1 = [qwait(kind='CMP')]
        seq1 = [qwait(kind='RAM', addr=0xff)]

    def test_compile(self):
        q1 = self.q1
        seq1 = [X90(q1), Y90(q1)]
        seq2 = [X(q1), Y(q1), Z(q1)]
        label(seq1)
        label(seq2)
        # mainLL, wfs1 = Compiler.compile_sequence(seq1 + seq2)
        seqIR = Compiler.compile_sequence(seq1 + seq2)
        # mainLL, wfs2 = Compiler.compile_sequence([X(q1), qif(0, seq1), Y(q1)])
        seqIR = Compiler.compile_sequence([X(q1), qif(0, seq1), Y(q1)])
        # assert(wfs1 == wfs2)
        # mainLL, wfs3 = Compiler.compile_sequence([X(q1), qif(0, seq1, seq2), Y(q1)])
        seqIR = Compiler.compile_sequence([X(q1), qif(0, seq1, seq2), Y(q1)])
        # assert(wfs1 == wfs3)


if __name__ == "__main__":

    unittest.main()
