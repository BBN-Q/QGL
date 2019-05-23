from QGL import *

def setup_test_lib():
    cl = ChannelLibrary(db_resource_name=":memory:")
    cl.clear()
    q1 = cl.new_qubit(label='q1')
    q2 = cl.new_qubit(label='q2')
    q3 = cl.new_qubit(label='q3')
    q4 = cl.new_qubit(label='q4')
    m1 = Measurement(label='M-q1', control_chan=q1, channel_db=cl.channelDatabase)
    m2 = Measurement(label='M-q2', control_chan=q2, channel_db=cl.channelDatabase)
    m3 = Measurement(label='M-q3', control_chan=q3, channel_db=cl.channelDatabase)
    m4 = Measurement(label='M-q4', control_chan=q4, channel_db=cl.channelDatabase)
    
    ChannelLibraries.channelLib.channelDict = {
        'q1': q1,
        'q2': q2,
        'q3': q3,
        'q4': q4,
        'q1q2': Edge(label='q1q2', source=q1, target=q2, channel_db=cl.channelDatabase),
        'q2q3': Edge(label='q2q3', source=q2, target=q3, channel_db=cl.channelDatabase),
        'q3q4': Edge(label='q3q4', source=q3, target=q4, channel_db=cl.channelDatabase)
    }
    cl.update_channelDict()
