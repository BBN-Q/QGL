'''
Copyright 2017 Raytheon BBN Technologies

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

from importlib import import_module
from pkgutil import walk_packages

import QGL.drivers

DRIVERS = dict()

for driver in walk_packages(QGL.drivers.__path__):
    DRIVERS[driver[1]] = import_module('QGL.drivers.' + driver[1])


def get_drivers():
    """
    Returns a name->module dictionary of device drivers known
    available
    """

    return DRIVERS


def add_driver(name, module):
    """
    Adds a new driver to the driver dictionary

    TODO: this may be incomplete; there will probably be more
    things we need to set up for each driver (for the plotter,
    if nothing else)
    """

    DRIVERS[name] = module
