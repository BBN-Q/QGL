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
from lark import Lark, tree, ast_utils, Transformer, v_args
from lark.tree import Meta
import ast

py="""
pots = 1000
pans = 2000
pens = 2001
fens = pens + 1
fens += 1
q1 = qubit()
"""
print("---- Python AST -----")
# print("Python AST:", tree)
print(ast.dump(ast.parse(py)))




this_module = sys.modules[__name__]

from dataclasses import dataclass

# Define AST passthroughs

class _Ast(ast_utils.Ast):
    pass

class _Decl(_Ast):
    pass

class _Expr(_Ast):
    pass

class _Assignment(_Ast):
    pass
# Define useful AST classes

@dataclass
class ConstDecl(_Decl):
    name: str
    value: _Expr

@dataclass
class Version(_Ast):
    number: int

@dataclass
class Assign(_Ast):
    value: object

@dataclass
class AugAssign(_Ast):
    operand: str
    value: object

class ClassicalAssignment(_Ast):
    def __init__(self, name, expr):
        self.name = name
        self.expr = expr
        print(name, expr)
        # return None
    def __str__(self):
        return(f"ClassicalAssignment(name={self.name}, expr={self.expr})")


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

    def run(self):
        transformer = ast_utils.create_transformer(sys.modules[__name__], QASM3Transformer())
        self.ast = transformer.transform(self.cst)
        return self.ast 
        # return QASM3Transformer().transform(self.cst)

    def ast_graph(self, filename):
        if _has_pydot:
            tree.pydot__tree_to_png(self.ast, filename)
        else:
            raise ModuleNotFoundError("Please install pydot to generate tree graphs.")


class QASM3Transformer(Transformer):
    def start(self, x):
        return x
    def expr(self, x):
        return x
    def NUMBER(self, tok):
        "Convert the value of `tok` from string to int, while maintaining line number & column."
        return float(tok) #tok.update(value=float(tok))
    def id(self, tok):
        return str(tok[0].value)


if __name__ == '__main__':
    q = QASM3Parser()
    with open(sys.argv[1]) as f:
        src = f.read()
        q.build_tree(src)
    q.cst_graph("qasm_tree.png")
    print("---- SRC -----")
    print(src)
    print("---- CST -----")
    print(q.cst)
    print("\n---- AST -----")
    # def prpr(x, n=1):
    #     if isinstance(x, list):

    for l in q.run():
        # if isinstance(l, list):
        print(l)
    # print(q.run())
    # q.ast_graph("qasm_tree_ast.png")
    
