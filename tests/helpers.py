from QGL import *

def setup_test_lib():
    q1 = Qubit(label='q1')
    q2 = Qubit(label='q2')
    q3 = Qubit(label='q3')
    q4 = Qubit(label='q4')

    ChannelLibrary.channelLib.channelDict = {
        'q1': q1,
        'q2': q2,
        'q3': q3,
        'q4': q4,
        'q1q2': Edge(label='q1q2', source=q1, target=q2),
        'q2q3': Edge(label='q2q3', source=q2, target=q3),
        'q3q4': Edge(label='q3q4', source=q3, target=q4)
    }
    ChannelLibrary.channelLib.build_connectivity_graph()
