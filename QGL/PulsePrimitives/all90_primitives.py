'''
Copyright 2016 Raytheon BBN Technologies

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''

from .common_primitives import overrideDefaults, _memoize, clear_pulse_cache, Id, Z, X90, X90m, Y90, Y90m, Z90, Z90m, Utheta, Xtheta, Ytheta, Ztheta, U90, U, arb_axis_drag,\
AC, DiAC, CNOT, CNOT_CR, flat_top_gaussian, echoCR, ZX90_CR, MEAS, MeasEcho, BLANK

#Setup the default 90/180 rotations
@_memoize
def X(qubit, **kwargs):
    return X90(qubit, **kwargs) + X90(qubit, **kwargs)

@_memoize
def Xm(qubit, **kwargs):
    return X90m(qubit, **kwargs) + X90m(qubit, **kwargs)

@_memoize
def Y(qubit, **kwargs):
    return Y90(qubit, **kwargs) + Y90(qubit, **kwargs)

@_memoize
def Ym(qubit, **kwargs):
    return Y90m(qubit, **kwargs) + Y90m(qubit, **kwargs)
