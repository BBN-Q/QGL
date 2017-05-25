import string


class BlockLabel(object):
    def __init__(self, label):
        self.label = label

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "{0}:".format(self.label)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.label == other.label
        return False

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.label)

    def promote(self, ptype):
        return self

    @property
    def length(self):
        return 0


def label(seq):
    if not isinstance(seq[0], BlockLabel):
        seq.insert(0, newlabel())
    return seq[0]


def endlabel(seq):
    if not isinstance(seq[-1], BlockLabel):
        seq.append(newlabel())
    return seq[-1]


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
