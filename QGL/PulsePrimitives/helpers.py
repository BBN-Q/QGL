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
from .. import PulseShapes
from .. import Channels
from .. import ChannelLibrary
import operator

from math import pi, sin, cos, acos, sqrt
import numpy as np
from ..PulseSequencer import Pulse, TAPulse, align
from functools import wraps, reduce

def overrideDefaults(chan, updateParams):
    '''Helper function to update any parameters passed in and fill in the defaults otherwise.'''
    # The default parameter list depends on the channel type so pull out of channel
    # Then update passed values
    paramDict = chan.pulseParams.copy()
    paramDict.update(updateParams)
    return paramDict


def _memoize(pulseFunc):
    """ Decorator for caching pulses to keep waveform memory usage down. """
    # namespacce the cache so we can access (and reset) from elsewhere
    _memoize.cache = {}

    @wraps(pulseFunc)
    def cacheWrap(*args, **kwargs):
        if kwargs:
            return pulseFunc(*args, **kwargs)
        key = (pulseFunc, args)
        if key not in _memoize.cache:
            _memoize.cache[key] = pulseFunc(*args)
        return _memoize.cache[key]

    return cacheWrap

def clear_pulse_cache():
    _memoize.cache = {}
