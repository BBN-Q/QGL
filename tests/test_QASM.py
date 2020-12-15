# Test QASM parsing and (eventually) compilation
import os
import unittest

from QGL.qasm.parse import QASM3Parser

def get_qasm_test_file(filename):
    path = os.path.dirname(os.path.abspath(__file__))
    file = os.path.join(path, "test_data", "qasm", filename)
    with open(file, "r") as f:
    	qasm = f.read()
    return qasm 

class ParseTestCase(unittest.TestCase):

	def setUp(self):
		pass

	def test_parse_simple(self):
		qasm = get_qasm_test_file("basic.qasm")
		parser = QASM3Parser()
		parser.build_tree(qasm)

if __name__ == "__main__":
	unittest.main()