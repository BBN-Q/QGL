from BlockLabel import label, endlabel

## QGL control-flow statements ##

def qif(condition, ifSeq, elseSeq):
	seqA = [GotoEq(condition, label(ifSeq))] + elseSeq
	seqB = ifSeq + [Goto(endlabel(elseSeq))]
	return seqA, seqB

def qwhile(condition, seq):
	return [GotoNeq(condition, endlabel(seq))] + seq

def qdowhile(condition, seq):
	return seq + qwhile(condition, seq)

## Sequencer primitives ##

class ControlInstruction(object):
	def __init__(self, instruction, condition=None, operator=None, target=None):
		self.instruction = instruction
		self.condition = condition
		self.operator = operator
		self.target = target

	def __repr__(self):
		return self.__str__()

	def __str__(self):
		result = self.instruction + "("
		if self.condition and self.operator:
			result += "{0} {1}, ".format(self.operator, self.condition)
		if self.target:
			result += str(self.target)
		result += ")"
		return result

def Goto(target):
	return ControlInstruction("GOTO", target=target)

def GotoEq(cond, target):
	return ControlInstruction("GOTO", cond, "==", target)

def GotoNeq(cond, target):
	return ControlInstruction("GOTO", cond, "!=", target)

def Call(target):
	return ControlInstruction("CALL", target=target)

def CallEq(cond, target):
	return ControlInstruction("CALL", cond, "==", target)

def CallNeq(cond, target):
	return ControlInstruction("CALL", cond, "!=", target)

def Return():
	return ControlInstruction("RETURN")

def ReturnEq(cond):
	return ControlInstruction("RETURN", cond, "==")

def ReturnNeq(cond):
	return ControlInstruction("RETURN", cond, "!=")
