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

import re, math
from sly import Lexer, Parser

class QASM3Lexer(Lexer):
	"""A lexer to tokenize the QASM3 language."""

    #set of QASM3 tokens
    tokens = {NUMCONST, STRCONST, BOOLCONST,
              VERSION, INCLUDE, MATH,
              IDENT, QUBIT, PHYSQUBIT, 
              PRECISION, WIDTH, SLICE, MATHFUNC,
              BIT, INTTYPE, UINTTYPE, FPTYPE, FLOATTYPE, ANGLETYPE, BOOLTYPE,
              CONSTTYPE ,LENGTHTYPE, STRETCHTYPE, TIME, LET,
              CNOT, GATE1, GPHASE, INV, POW, CTRL, RESET, MEAS,
              BITOP, BOOLOP, NUMOP, ASSIGN,
              IF, ELSE, FOR, WHILE, CONTINUE, BREAK, END,
              KERNEL, DEF, PRAGMA, LENOF, BOX, BARRIER} 
    
    #Ignored characters
    ignore = ' \t'
    
    #Ignore comments
    ignore_single_comment = r'/{2}.*'
    ignore_multi_comment  = r'/\*[\s\S]*?\*/'
    
    @_(r'\n+')
    def ignore_newline(self, t):
        self.lineno += t.value.count('\n')
    
    literals = {':', ';', ',', '=', '(', ')', '{', '}', '@'}
    
    @_(r'true|false')
    def BOOLCONST(self, t):
    	"""Matches a boolean constant true/false."""
        t.value = True if t.value == "true" else False
        return t
    
    @_(r'pi', r'tau', r'e', r'[+-]?0b[01]+', r'[+-]?0x[0-9a-fA-F]+', 
    	r'[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?')
    def NUMCONST(self, t):
        """Matches a numeric constant as either binary (0b), hex (0x), or generic floating point."""
        if t.value == 'pi':
            t.value = math.pi
        elif t.value == 'tau':
            t.value = 2.0*math.pi
        elif t.value == 'e':
            t.value = math.e
        else:
            if t.value.startswith('0x'):
                t.value = int(t.value[2:], 16)
            elif t.value.startswith('0b'):
                t.value= int(t.value[2:], 2)
            t.value = float(t.value)
        return t
    
    @_(r'\"[\S\t\f ]*\"')
    def STRCONST(self, t):
       """Matches a string literal constant."""
       t.value = str(t.value[1:-1])
       return t

    VERSION = r'QASMVERSION'		
    INCLUDE = r'include'			#Include another file
    
    QUBIT     = r'qubit|qreg'		#Qubit type
    PHYSQUBIT = r'\%\d+'            #Physical qubit 
    BIT       = r'bit|creg'         #Classical bit/register identifier
    
    PRECISION = r'\d+:\d+:\d+'      #Precision identifier for fixed-point numbers
    
    INTTYPE    = r'int'				
    UINTTYPE   = r'uint'
    FPTYPE     = r'fixed'
    FLOATTYPE  = r'float'
    ANGLETYPE  = r'angle'
    BOOLTYPE   = r'bool'
    CONSTTYPE  = r'const'
    LENGTHTYPE = r'length'
    
    TIME = r'dt|ns|us|ms|s'
    
    LET = r'let'
    
    CNOT   = r'CX'
    GATE1  = r'U'
    GPHASE = r'gphase'
    INV    = r'inv'
    POW    = r'pow'
    CTRL   = r'ctrl'
    
    RESET  = r'reset'
    MEAS   = r'measure'
    
    ASSIGN = r'->'
    
    #These are the operations on bitstring, booleans, and numeric types
    BITOP  = r'&|\||\^|<<|>>|~|popcount|rotl|rotr'
    BOOLOP = r'[><]=?|==|!=?|&&|\|\||in'
    NUMOP  = r'\+[\+=]?|-[-=]?|\*=|/=?'
        
    IF       = r'if'
    ELSE     = r'else'
    FOR      = r'for'
    WHILE    = r'while'
    CONTINUE = r'continue'
    BREAK    = r'break'
    END      = r'end'
    KERNEL   = r'kernel'
    DEF      = r'def'
    LENOF    = r'lengthof'
    BARRIER  = r'barrier'
    
    @_(r'boxas|boxto')
    def BOX(self, t):
        t.value = t.value[-2:]
        return t
    
    @_(r'stretch\d{0,3}')
    def STRETCHTYPE(self, t):
        match = re.search(r'\d+', t.value)
        if match:
            t.value = int(match.group())
        else:
            t.value = 0
        return t
    
    #Built-in math functions
    @_(r'sqrt|floor|ceiling|log|pow|div|mod|sin|cos|tan')
    def MATHFUNC(self, t):
        if t.value == "ceiling":
            t.value = math.ceil
        else:
            t.value = getattr(math, t.value)
        return t
        
    @_(r'\[\d+\]')
    def WIDTH(self, t):
       t.value = int(t.value[1:-1])
       return t
    
    @_(r'\[\d+:\d+\]|\[\d+:\d+:\d+\]')
    def SLICE(self, t):
        match = re.match(r'\[(\d+):(\d+)\]|\[(\d+):(\d+):(\d+)\]', t.value)
        t.value = [int(g) for g in match.groups() if g]
        return t
    
    @_(r'\#PRAGMA[\S\t\f ]+')
    def PRAMGA(self, t):
        match = re.match(r'\#PRAGMA([\S\t\f ]+)', t.value)
        t.value = match.groups()[0].lstrip()
        return t
    
    IDENT = r'[a-zA-Z_%][a-zA-Z0-9_]*' #Variable identifier
    
    def error(self, t):
        print('Line %d: Bad character %r' % (self.lineno, t.value[0]))
        self.index += 1
