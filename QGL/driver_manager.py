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

DRIVERS = dict()
LOADED = False

DEFAULT_DEVICES = [
        'APS2Pattern', 'APSPattern', 'TekPattern'
        ]

def get_drivers():
    """
    Returns a name->module dictionary of device drivers known
    available
    """

    global LOADED

    if not LOADED:
        for dev_name in DEFAULT_DEVICES:
            DRIVERS[dev_name] = import_module('QGL.drivers.' + dev_name)
        LOADED = True

    return DRIVERS


def add_driver(name, module):
    """
    Adds a new driver to the driver dictionary

    TODO: check whether there is already a device loaded with
    the given name, and at least warn the user
    """

    DRIVERS[name] = module
