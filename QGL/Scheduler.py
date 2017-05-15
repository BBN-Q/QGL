
from .Channels import Edge
from .PulseSequencer import PulseBlock
from .ControlFlow import Barrier, ControlInstruction
from .PatternUtils import flatten
from warnings import warn

def schedule(seq):
    '''
    Takes in a "serial" sequence of operations and creates a new sequence which
    packs those operations in a maximally concurrent fashion. The user inserts
    `Barrier()` operations to prevent the scheduler from moving operations "earlier"
    than the barrier.
    '''

    # dictionary of per-qubit time counters
    counters = {}
    out_seq = []

    for instr in flatten(seq):
        if isinstance(instr, Barrier):
            synchronize_counters(counters, instr.chanlist)
            continue
        channels = get_channels(instr)
        # find the most advanced counter in the channel set
        idx = max(counters.get(ch, 0) for ch in channels)

        if idx > len(out_seq) - 1:
            out_seq.append(instr)
        else:
            out_seq[idx] *= instr
        # advance the channel counter(s)
        for ch in channels:
            counters[ch] = idx + 1

    return out_seq

def get_channels(instr):
    '''
    Normalizes the various ways 'channels' can be encoded in instruction into
    a tuple of channels
    '''
    if not hasattr(instr, 'channel'):
        warn("instruction %s does not have a 'channel' property", instr)
        return None
    if isinstance(instr, PulseBlock):
        return tuple(instr.channel)
    elif isinstance(instr.channel, Edge):
        return (instr.channel.source, instr.channel.target)
    else:
        return (instr.channel,)

def synchronize_counters(counters, channels):
    '''
    Advance the counter for each channel in 'channels' to the largest count
    in the set.
    '''
    max_idx = max(counters.get(ch, 0) for ch in channels)
    for ch in channels:
        counters[ch] = max_idx
