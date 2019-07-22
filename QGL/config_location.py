# -*- coding: UTF-8 -*-
'''
Created on Jul 24, 2017

@author: dan.ellard@raytheon.com

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

HOW TO USE:

The module or file that is serving as the "main" should begin
with something like the following, BEFORE the importing of
any of the QGL files:

from . import config_location
config_location.config(PATH_TO_CONFIG_FILE)
from . inport config

If you don't wish to set the path to the configuration file
explicitly, then the first two lines should be omitted.
The last line may be unnecessary because many of the files in
this directory begin with something like "from QGL import *",
which imports config.

NOTE: This file must be imported, and any defaults overridden,
before config.py is imported.  When config.py is imported, it
will use the get_config_path method from this module to find the
location of the JSON config file to use.  After this happens,
the config variables are set and will not be altered by changing
the CONFIG_PATH or the _CONFIG_FILE_NAME (and they probably
SHOULD NOT be modified, because some of the values taken from
the config file directly change the behavior of the other modules,
so changing things on fly may result in an inconsistent state).

'''

import os.path

CONFIG_PATH = None

_CONFIG_FILE_NAME = 'BBN_config.py'

def _set_default_config_path():
    """
    Set the path for the "default" configuration file:

    1. if there's one in the current working directory, then use it

    2. otherwise, if there's one in dirname(__file__), then use it

    3. otherwise, leave the path at None and assume that the user will
        set it explicitly instead of relying on a default
    """

    cwd_path = os.path.join(os.getcwd(), _CONFIG_FILE_NAME)
    def_path = os.path.join(os.path.dirname(__file__), _CONFIG_FILE_NAME)
    
    # We could also check whether the file is readable, etc.
    # I think the error messages make more sense if we let this
    # fail later, but I don't know what the users will think is
    # the most intelligible.

    if os.path.isfile(cwd_path):
        return cwd_path
    elif os.path.isfile(def_path):
        return def_path
    else:
        return None


def config(path):
    """
    Set CONFIG_FILE, the path to the configuration JSON file,
    to the given path, and return the new path.
    
    No error or sanity checking is performed.
    """

    global CONFIG_FILE

    CONFIG_FILE = path
    return CONFIG_FILE


def set_config_filename(filename):
    """
    Change the "default" filename for the config file from the
    current value of CONFIG_FILE_NAME to the given filename string.

    The filename is assumed to be a simple filename (no directory
    components).  Partial directory paths might work, or might not.
    No sanity checking is done on the filename.

    This function has no useful effect if CONFIG_FILE is non-None;
    it is used only to compute the CONFIG_FILE, and if the CONFIG_FILE
    has already been set, then CONFIG_FILE will not change as a
    result of changing the CONFIG_FILE_NAME.
    """

    global _CONFIG_FILE_NAME

    _CONFIG_FILE_NAME = filename
    return filename


def get_config_path():
    """
    Retrieve the current value of the path to the configuration file
    """

    global CONFIG_PATH

    if not CONFIG_PATH:
        CONFIG_PATH = _set_default_config_path()

    return CONFIG_PATH
