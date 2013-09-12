import string

class BlockLabel(object):
	def __init__(self, label, offset=0):
		self.label = label
		self.offset = offset

	def __add__(self, offset):
		# allows for addition/subtraction of offsets to labeled blocks
		return BlockLabel(self.label, self.offset + offset)

	def __sub__(self, offset):
		return self.__add__(-offset)

	def __repr__(self):
		return self.__str__()

	def __str__(self):
		if self.offset == 0:
			return self.label
		else:
			return "{0}+{1}".format(self.label, self.offset)

	@property
	def empty(self):
		if self.label == "" and self.offset == 0:
			return True
		else:
			return False

def label(seq):
	# label's are attached to the first block in a sequence
	# make sure we have a PulseBlock
	seq[0] = seq[0].promote()
	if seq[0].label == None or seq[0].label.empty:
		seq[0].label = newlabel()
	return seq[0].label

def endlabel(seq):
	return label(seq) + len(seq)

numlabels = 0
def newlabel():
	labelSet = string.ascii_uppercase
	global numlabels
	numlabels += 1
	if numlabels > len(labelSet):
		NameError("Ran out of labels")
	return BlockLabel("seq" + labelSet[numlabels-1])
