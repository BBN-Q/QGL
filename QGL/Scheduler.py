
from .Channels import Edge, Measurement
from .PulseSequencer import PulseBlock
from .ControlFlow import Barrier, ControlInstruction
from .BlockLabel import BlockLabel
from .PatternUtils import flatten
from .ChannelLibraries import QubitFactory
from warnings import warn

def schedule(seq):
    '''
    Takes in a "serial" sequence of operations and creates a new sequence which
    packs those operations in a maximally concurrent fashion. The user inserts
    `Barrier()` operations to prevent the scheduler from moving operations "earlier"
    than the barrier.
    '''

    # dictionary of per-qubit time counters
    # each value will encode the last filled time bin per channnel (key)
    counters = {}

    out_seq = []
    channel_set = find_all_channels(seq)

    for instr in seq:
        if isinstance(instr, Barrier):
            synchronize_counters(counters, instr.chanlist)
            continue

        channels = get_channels(instr, channel_set)
        # find the most advanced counter in the channel set
        idx = max(counters.get(ch, 0) for ch in channels)

        if (idx >= len(out_seq)) or isinstance(out_seq[idx], ControlInstruction):
            out_seq.append(instr)
        else:
            out_seq[idx] *= instr
        # advance the channel counter(s)
        for ch in channels:
            counters[ch] = idx + 1

    return out_seq

def get_channels(instr, channel_set=None):
    '''
    Normalizes the various ways 'channels' can be encoded in instructions into
    a tuple of channels
    TODO this suggests that we should normalize the names of these properties
    '''
    if isinstance(instr, PulseBlock):
        return tuple(instr.channel)
    elif isinstance(instr, (ControlInstruction, BlockLabel)):
        # these instruction types are assumed to broadcast
        return channel_set
    elif isinstance(instr, Barrier):
        return chanlist
    elif not hasattr(instr, 'channel'):
        warn("instruction %s does not have a 'channel' property", instr)
        return None
    elif isinstance(instr.channel, Edge):
        return (instr.channel.source, instr.channel.target)
    elif isinstance(instr.channel, Measurement):
        # TODO update Measurement channels to contain a reference back to their
        # parent Qubit
        _, qubit_name = instr.channel.label.rsplit('-', 1)
        return (QubitFactory(qubit_name),)
    else:
        return (instr.channel,)

def find_all_channels(seq):
    channels = set()
    for instr in seq:
        instr_channels = get_channels(instr)
        if instr_channels is None:
            continue
        for ch in instr_channels:
            channels.add(ch)
    return channels

def synchronize_counters(counters, channels):
    '''
    Advance the counter for each channel in 'channels' to the largest count
    in the set.
    '''
    max_idx = max(counters.get(ch, 0) for ch in channels)
    for ch in channels:
        counters[ch] = max_idx
