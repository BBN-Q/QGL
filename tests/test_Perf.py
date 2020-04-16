#simple file added for performace testing...

from QGL import *
import random

def random_sequence_test():
    N = 1000

    cl = ChannelLibrary()
    q1 = QubitFactory("q1")

    gates = (X90(q1), Y90(q1), X(q1), Y(q1), Id(q1, length=1e-7))

    seqs = []
    ks = random.choices(range(1,100), k = N)
    for k in ks:
        seqs.append(random.choices(gates, k=k) + [MEAS(q1)])

    compile_to_hardware(seqs, 'test/test')
