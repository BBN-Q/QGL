from BlockLabel import label, endlabel
from functools import wraps


## QGL control-flow statements ##

def qif(mask, ifSeq, elseSeq):
	endlabel(elseSeq) # make sure to populate label of elseSeq before using it
	seqA = [CmpEq(mask), Goto(label(ifSeq))] + elseSeq
	seqB = ifSeq + [Goto(endlabel(elseSeq))]
	return seqA, seqB

def qwhile(mask, seq):
	return [CmpNeq(mask), Goto(endlabel(seq))] + seq

def qdowhile(mask, seq):
	return seq + [CmpEq(mask), Goto(label(seq))]

def _qfunction(func):
	@wraps(func)
	def crfunc(*args):
		if not crfunc.target:
			crfunc.seq = func(*args)
			crfunc.target = label(crfunc.seq)
		return [Call(crfunc.target)], crfunc.seq + [Return()]
	crfunc.target = None
	crfunc.seq = None
	return crfunc

## Sequencer primitives ##

class ComparisonInstruction(object):
	def __init__(self, mask, operator):
		self.mask = mask
		self.operator = operator

	def __repr__(self):
		return self.__str__()

	def __str__(self):
		return "CMP " + self.operator + " " + str(self.mask)

def CmpEq(mask):
	return ComparisonInstruction(mask, "==")

def CmpNeq(mask):
	return ComparisonInstruction(mask, "!=")

def CmpLt(mask):
	return ComparisonInstruction(mask, "<")

def CmpGt(mask):
	return ComparisonInstruction(mask, ">")

class ControlInstruction(object):
	def __init__(self, instruction, target=None):
		self.instruction = instruction
		self.target = target

	def __repr__(self):
		return self.__str__()

	def __str__(self):
		result = self.instruction
		if self.target:
			result += "(" + str(self.target) + ")"
		return result

def Goto(target):
	return ControlInstruction("GOTO", target=target)

def Call(target):
	return ControlInstruction("CALL", target=target)

def Return():
	return ControlInstruction("RETURN")
