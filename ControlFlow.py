from BlockLabel import newlabel, label, endlabel
from PulseSequencer import Pulse
from functools import wraps
from mm import multimethod

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
	return qfunction_seq[target]

@multimethod(int, Pulse)
def repeat(n, p):
    p.repeat = round(n)
    return p

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
		elif self.value:
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
	def __init__(self, target):
		super(Goto, self).__init__("GOTO", target=target)

class Call(ControlInstruction):
	def __init__(self, target):
		super(Call, self).__init__("CALL", target=target)

class Return(ControlInstruction):
	def __init__(self):
		super(Return, self).__init__("RETURN")

class LoadRepeat(ControlInstruction):
	def __init__(self, n):
		super(LoadRepeat, self).__init__("LOAD", value=n)

class Repeat(ControlInstruction):
	def __init__(self, target):
		super(Repeat, self).__init__("REPEAT", target=target)

class Wait(ControlInstruction):
	def __init__(self):
		super(Wait, self).__init__("WAIT")

class LoadCmp(ControlInstruction):
	def __init__(self):
		super(LoadCmp, self).__init__("LOADCMP")

class Sync(ControlInstruction):
	def __init__(self):
		super(Sync, self).__init__("SYNC")

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
