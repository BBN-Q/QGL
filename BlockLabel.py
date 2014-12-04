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

	def __eq__(self, other):
		if isinstance(other, self.__class__):
			return self.__dict__ == other.__dict__
		return False

	def __hash__(self):
		return hash(self.label) ^ hash(self.offset)

def label(seq):
	# label's are attached to the first block in a sequence
	# make sure we have a PulseBlock
	seq[0] = seq[0].promote()
	if seq[0].label == None:
		seq[0].label = newlabel()
	return seq[0].label

def endlabel(seq):
	return label(seq) + len(seq)

def newlabel():
	label = BlockLabel(asciibase(newlabel.numlabels))
	newlabel.numlabels += 1
	return label

newlabel.numlabels = 0

def asciibase(x):
	'''
	Converts x to base 26 composed of the uppercase alphabet
	'''
	digits = string.ascii_uppercase
	s = ''
	while 1:
		s = digits[x % 26] + s
		x = x // 26
		if x == 0:
			break
	return s
