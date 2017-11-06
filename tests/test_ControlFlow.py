import unittest

from QGL import *
from QGL.ControlFlow import CmpEq, CmpNeq, Goto, Call, Return, LoadRepeat, Repeat
from QGL.BlockLabel import label, endlabel


class ControlFlowTest(unittest.TestCase):
    def setUp(self):
        self.q1 = Qubit(label='q1')
        self.q2 = Qubit(label='q2')
        
        ChannelLibrary(library_file=None) # Create a blank ChannelLibrary
        ChannelLibraries.channelLib.channelDict = {'q1': self.q1, 'q2': self.q2}
        ChannelLibraries.channelLib.build_connectivity_graph()

    def test_qif(self):
        q1 = self.q1
        seq1 = [X90(q1), Y90(q1)]
        seq2 = [X(q1), Y(q1), Z(q1)]
        label(seq1)
        label(seq2)
        # print qif(0, seq1, seq2)
        # print ([CmpEq(0), Goto(label(seq1))] + seq2 + [Goto(endlabel(seq1))] + seq1
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
        seq1 = [qwait(), qwait(kind="CMP")]
        assert (isinstance(seq1[0], ControlFlow.Wait))
        assert (isinstance(seq1[1], ControlFlow.LoadCmp))

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
