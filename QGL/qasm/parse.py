"""
Original Author: Guilhem Ribeill

Copyright 2020 Raytheon BBN Technologies

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
import os, sys
import logging 
from lark import Lark, tree, logger 

#logger.setLevel(logging.DEBUG)

try:
    import pydot
    _has_pydot = True
except ImportError:
    _has_pydot = False

grammar_path = os.path.join(os.path.dirname(
                            os.path.abspath(__file__)),
                             "grammar.lark")
with open(grammar_path, "r") as f:
    _QASM_GRAMMAR = f.read()

class QASM3Parser(Lark):

    def __init__(self, **kwargs):
        super().__init__(_QASM_GRAMMAR, **kwargs)
        self.cst = None

    def build_tree(self, input):
        self.cst = self.parse(input)

    def __str__(self):
        return self.cst.pretty()

    def cst_graph(self, filename):
        if _has_pydot:
            tree.pydot__tree_to_png(self.cst, filename)
        else:
            raise ModuleNotFoundError("Please install pydot to generate tree graphs.")

if __name__ == "__main__":
    with open(sys.argv[1]) as f:
        parser = QASM3Parser(debug=True)
        parser.build_tree(f.read())
        print(parser)
        parser.cst_graph("test_qasm.png")
