from BlockLabel import label, endlabel
from functools import wraps


## QGL control-flow statements ##

def qif(mask, ifSeq, elseSeq=None):
	if elseSeq:
		endlabel(elseSeq) # make sure to populate label of elseSeq before using it
		return [CmpEq(mask), Goto(label(ifSeq))] + [(elseSeq, ifSeq + [Goto(endlabel(elseSeq))])]
	else:
		endlabel(ifSeq)
		return [CmpNeq(mask), Goto(endlabel(ifSeq))] + ifSeq

def qwhile(mask, seq):
	return [CmpNeq(mask), Goto(endlabel(seq))] + seq

def qdowhile(mask, seq):
	return seq + [CmpEq(mask), Goto(label(seq))]

def qfunction(func):
	# caches for sequences and labels
	seq = {}
	target = {}
	@wraps(func)
	def crfunc(*args):
		if args not in target:
			seq[args] = func(*args)
			target[args] = label(seq[args])
		return [Call(target[args])], seq[args] + [Return()]
	return crfunc

## Sequencer primitives ##

class ComparisonInstruction(object):
	def __init__(self, mask, operator):
		self.mask = mask
		self.operator = operator
		self.label = None

	def __repr__(self):
		return self.__str__()

	def __str__(self):
		labelPart = "{0}: ".format(self.label) if self.label else ""
		return labelPart + "CMP " + self.operator + " " + str(self.mask)

	def __eq__(self, other):
		return self.__dict__ == other.__dict__

	def promote(self):
		return self

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
		self.label = None

	def __repr__(self):
		return self.__str__()

	def __str__(self):
		labelPart = "{0}: ".format(self.label) if self.label else ""
		result = labelPart + self.instruction
		if self.target:
			result += "(" + str(self.target) + ")"
		return result

	def __eq__(self, other):
		return self.__dict__ == other.__dict__

	def promote(self):
		return self

def Goto(target):
	return ControlInstruction("GOTO", target=target)

def Call(target):
	return ControlInstruction("CALL", target=target)

def Return():
	return ControlInstruction("RETURN")
