from .BlockLabel import newlabel, label, endlabel
from .PulseSequencer import Pulse
from functools import wraps
from .mm import multimethod

## QGL control-flow statements ##

# FIXME qif, qwhile, and qdowhile need an address to pass to CMP

def qif(value, ifSeq, elseSeq=None):
    '''
    Create a branching sequence that
    '''
    if elseSeq:
        return [CmpEq("m", value), Goto(label(ifSeq))] + elseSeq + [
            Goto(endlabel(ifSeq))
        ] + ifSeq
    else:
        endlabel(ifSeq)
        return [CmpNeq("m", value), Goto(endlabel(ifSeq))] + ifSeq


def qwhile(value, seq):
    label1 = newlabel()
    label2 = newlabel()
    return [label1, CmpNeq("m", value), Goto(label2)] + seq + [Goto(label1), label2]


def qdowhile(value, seq):
    return seq + [CmpEq("m", value), Goto(label(seq))]

# caches for sequences and labels
qfunction_seq = {}


def qfunction(func):
    target = {}

    @wraps(func)
    def crfunc(*args, **kwargs):
        if args not in target:
            seq = func(*args, **kwargs) + [Return()]
            target[args] = label(seq)
            qfunction_seq[label(seq)] = seq
        return Call(target[args])

    return crfunc


def qfunction_specialization(target):
    # QGL2 uses Call with a label, pointing elsewhere in the same
    # sequenece, wanting control passed there.
    # In those cases, there was no qfunction registering
    # a sequence that needs to be inserted; the pulses
    # are already in a sequence, so are not in this hash.
    if target in qfunction_seq:
        return qfunction_seq[target]
    else:
        # QGL2 case
        return list()


@multimethod(int, Pulse)
def repeat(n, p):
    return repeat(n, [p])


@multimethod(int, list)
def repeat(n, seq):
    if n < 1:
        return None
    elif n == 1:
        return seq
    else:
        label(seq)
        return [LoadRepeat(n)] + seq + [Repeat(label(seq))]


# utility to repeat all sequences the same number of times
def repeatall(n, seqs):
    for ct in range(len(seqs)):
        seqs[ct] = repeat(n, seqs[ct])
    return seqs


def qwait(channels=None, kind="TRIG"):
    '''
    Insert a WAIT or LOADCMP command on the target channels (an iterable, or None)
    '''
    if kind == "TRIG":
        return Wait(channels)
    else:
        return LoadCmp(channels)


def qsync(channels=None):
    '''
    Insert a SYNC command on the target channels (an iterable, or None).
    '''
    return Sync(channels)

## Sequencer primitives ##

class ControlInstruction(object):
    def __init__(self, instruction, channels=None, target=None, value=None):
        self.channels = channels
        self.instruction = instruction
        self.target = target  #refactor into payload field??
        self.value = value
        self.length = 0

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        result = self.instruction + "("
        chan_str = str(self.channels) if self.channels else None
        target_str = str(self.target) if self.target else None
        value_str = str(self.value) if self.value else None
        result += ", ".join(filter(None, [chan_str, target_str, value_str]))
        result += ")"
        return result

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        return not self == other

    def promote(self, ptype):
        return self

class Store(ControlInstruction):
    def __init__(self, target, value):
        super(Store, self).__init__("STORE", target=target, value=value)

class Goto(ControlInstruction):
    # target is a BlockLabel
    def __init__(self, target):
        super(Goto, self).__init__("GOTO", target=target)


class Call(ControlInstruction):
    # target is a BlockLabel
    def __init__(self, target):
        super(Call, self).__init__("CALL", target=target)


class Return(ControlInstruction):
    def __init__(self):
        super(Return, self).__init__("RETURN")


class LoadRepeat(ControlInstruction):
    # n is an integer
    def __init__(self, n):
        super(LoadRepeat, self).__init__("LOAD", value=n)


class Repeat(ControlInstruction):
    # target is a BlockLabel
    def __init__(self, target):
        super(Repeat, self).__init__("REPEAT", target=target)


class Wait(ControlInstruction):
    def __init__(self, channels=None):
        super(Wait, self).__init__("WAIT", channels)


class LoadCmp(ControlInstruction):
    def __init__(self, channels=None):
        super(LoadCmp, self).__init__("LOADCMP", channels)


class Sync(ControlInstruction):
    def __init__(self, channels=None):
        super(Sync, self).__init__("SYNC", channels)


class ComparisonInstruction(ControlInstruction):
    def __init__(self, operator, address, value):
        super(ComparisonInstruction, self).__init__("CMP")
        self.operator = operator
        self.address = address
        self.value = value

    def __str__(self):
        return "CMP " + str(self.address) + " " + self.operator + " " + \
            str(self.value)


def CmpEq(address, value):
    return ComparisonInstruction("==", address, value)


def CmpNeq(address, value):
    return ComparisonInstruction("!=", address, value)


def CmpLt(address, value):
    return ComparisonInstruction("<", address, value)


def CmpGt(address, value):
    return ComparisonInstruction(">", address, value)

class Barrier(ControlInstruction):
    def __init__(self, *chanlist):
        # chanlist is a list of Channel instances
        super(Barrier, self).__init__("BARRIER")
        self.chanlist = chanlist

    def __str__(self):
        base = "BARRIER({0})".format(self.chanlist)
        return base
