from .BlockLabel import newlabel, label, endlabel
from .PulseSequencer import Pulse
from functools import wraps
from .mm import multimethod

## QGL control-flow statements ##

def qif(mask, ifSeq, elseSeq=None):
	if elseSeq:
		return [CmpEq(mask), Goto(label(ifSeq))] + elseSeq + [Goto(endlabel(ifSeq))] + ifSeq
	else:
		endlabel(ifSeq)
		return [CmpNeq(mask), Goto(endlabel(ifSeq))] + ifSeq

def qwhile(mask, seq):
	label1 = newlabel()
	label2 = newlabel()
	return [label1, CmpNeq(mask), Goto(label2)] + seq + [Goto(label1), label2]

def qdowhile(mask, seq):
	return seq + [CmpEq(mask), Goto(label(seq))]

# caches for sequences and labels
qfunction_seq = {}
def qfunction(func):
	target = {}
	@wraps(func)
	def crfunc(*args):
		if args not in target:
			seq = func(*args) + [Return()]
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

def qwait(kind="TRIG"):
	if kind == "TRIG":
		return Wait()
	else:
		return LoadCmp()

def qsync():
	return Sync()

## Sequencer primitives ##

class ControlInstruction(object):
	def __init__(self, instruction, target=None, value=None):
		self.instruction = instruction
		self.target = target #refactor into payload field??
		self.value = value
		self.label = None
		self.length = 0

	def __repr__(self):
		return self.__str__()

	def __str__(self):
		labelPart = "{0}: ".format(self.label) if self.label else ""
		result = labelPart + self.instruction
		if self.target:
			result += "(" + str(self.target) + ")"
		elif self.value or self.value == 0:
			result += "(" + str(self.value) + ")"
		return result

	def __eq__(self, other):
		if isinstance(other, self.__class__):
			# ignore label in equality testing
			mydict = self.__dict__.copy()
			otherdict = other.__dict__.copy()
			mydict.pop('label')
			otherdict.pop('label')
			return mydict == otherdict
		return False

	def __ne__(self, other):
		return not self == other

	def promote(self):
		return self

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
	def __init__(self):
		super(Wait, self).__init__("WAIT")

# FIXME: This is not supported by the hardware yet
# This is a Wait that should only wait for the channels
# listed in chanlist (not all channels, as in Wait)
class WaitSome(ControlInstruction):
        # chanlist is a list of Channel instances
	def __init__(self, chanlist):
                # Until HW really supports waitsome, use wait so things compile better
		super(WaitSome, self).__init__("WAIT")
#		super(WaitSome, self).__init__("WAITSOME")
		# The channels to wait on
		self.chanlist = chanlist

	def __str__(self):
		base = super(WaitSome, self).__str__()
		base += " on Channels: %s" % str(self.chanlist)
		return base

class LoadCmp(ControlInstruction):
	def __init__(self):
		super(LoadCmp, self).__init__("LOADCMP")

class Sync(ControlInstruction):
	def __init__(self):
		super(Sync, self).__init__("SYNC")

# For use within QGL2 only
# Marks a barrier at start or end of blocks
# that can be run concurrently (with concur blocks)
# Should be replaced by QGL2 with the proper length Id pulse,
# or a Sync then a Wait if the block is of indeterminate length.
class Barrier(ControlInstruction):
        # chanlist is a list of Channel instances
        # ctr is an opaque string, unique per channel
        # (but appearing once for the program for each channel
        # in chanlist)
        def __init__(self, ctr, chanlist):
                super(Barrier, self).__init__("BARRIER", value=ctr)
                # Consider adding a start/end marker,
                # to help line up barriers across sequences.
                # FIXME: self.start
                self.chanlist = chanlist
        def __str__(self):
                base = super(Barrier, self).__str__()
                base += " on Channels: %s" % str(self.chanlist)
                return base

class ComparisonInstruction(ControlInstruction):
	def __init__(self, mask, operator):
		super(ComparisonInstruction, self).__init__("CMP")
		self.mask = mask
		self.operator = operator

	def __str__(self):
		labelPart = "{0}: ".format(self.label) if self.label else ""
		return labelPart + "CMP " + self.operator + " " + str(self.mask)

def CmpEq(mask):
	return ComparisonInstruction(mask, "==")

def CmpNeq(mask):
	return ComparisonInstruction(mask, "!=")

def CmpLt(mask):
	return ComparisonInstruction(mask, "<")

def CmpGt(mask):
	return ComparisonInstruction(mask, ">")
