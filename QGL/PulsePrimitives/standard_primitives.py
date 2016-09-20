'''
Copyright 2013 Raytheon BBN Technologies

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

from .helpers import overrideDefaults, _memoize, clear_pulse_cache
from .common_primitives import Xtheta, Ytheta
#Setup the default 90/180 rotations
@_memoize
def X(qubit, **kwargs):
    return Xtheta(qubit,
                  qubit.pulseParams['piAmp'],
                  label="X",
                  ignoredStrParams=['amp'],
                  **kwargs)

@_memoize
def Xm(qubit, **kwargs):
    return Xtheta(qubit,
                  -qubit.pulseParams['piAmp'],
                  label="Xm",
                  ignoredStrParams=['amp'],
                  **kwargs)

@_memoize
def Y(qubit, **kwargs):
    return Ytheta(qubit,
                  qubit.pulseParams['piAmp'],
                  label="Y",
                  ignoredStrParams=['amp'],
                  **kwargs)

@_memoize
def Ym(qubit, **kwargs):
    return Ytheta(qubit,
                  -qubit.pulseParams['piAmp'],
                  label="Ym",
                  ignoredStrParams=['amp'],
                  **kwargs)
