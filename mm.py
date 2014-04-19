'''
Taken from Guido van Rossum's Five-Minute MultiMethods in Python
http://www.artima.com/weblogs/viewpost.jsp?thread=101605

Examples:
from mm import multimethod

@multimethod(int, int)
def foo(a, b):
    ...code for two ints...

@multimethod(float, float):
def foo(a, b):
    ...code for two floats...

@multimethod(str, str):
def foo(a, b):
    ...code for two strings...
'''

registry = {}

class MultiMethod(object):
    def __init__(self, name):
        self.name = name
        self.typemap = {}
    def __call__(self, *args):
        types = tuple(arg.__class__ for arg in args) # a generator expression!
        function = self.typemap.get(types)
        if function is None:
            raise TypeError("no match")
        return function(*args)
    def register(self, types, function):
        if types in self.typemap:
            raise TypeError("duplicate registration")
        self.typemap[types] = function

def multimethod(*types):
    def register(function):
        name = function.__name__
        mm = registry.get(name)
        if mm is None:
            mm = registry[name] = MultiMethod(name)
        mm.register(types, function)
        return mm
    return register
