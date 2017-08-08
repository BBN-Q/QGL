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

The rule for choosing the config file is if a config file path
is specified explicitly (using the method shown below) then
that path is used; otherwise, if the QGLCFGFILE environment
variable is set, its value is used as the path; otherwise
the current working directory is checked, and finally the
"default" directory where the module resides is checked.

To set the path explicitly, the module or file that is serving
as the "main" should begin with something like the following,
BEFORE the importing of any of the other QGL files:

# The import path may depend on other things and require
# additional qualification
import qgl_config_loc
qgl_config_loc.config(PATH_TO_CONFIG_FILE)

# and then the rest of the imports, as usual

If you don't wish to set the path to the configuration file
explicitly, then these lines should be omitted.

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

import os

CONFIG_PATH = None

_CONFIG_FILE_NAME = 'config.json'

def _set_default_config_path():
    """
    Set the path for the "default" configuration file:

    1. if there's a QGLCFGFILE environment variable set, return it

    2. if there's a _CONFIG_FILE_NAME in the current working directory,
        then return the path to it

    3. otherwise, return dirname(__file__) + '/' + _CONFIG_FILE_NAME

    Note that this routine does no error checking about
    whether or not the return value is usable.  If the file doesn't
    exist (in cases 1 or 3) it will be created from a template later.
    """

    env_path = os.getenv('QGLCFGFILE')
    cwd_path = os.path.join(os.getcwd(), _CONFIG_FILE_NAME)
    def_path = os.path.join(
            os.path.dirname(__file__), 'QGL', _CONFIG_FILE_NAME)

    # We could also check whether the file is readable, etc.
    # I think the error messages make more sense if we let this
    # fail later, but I don't know what the users will think is
    # the most intelligible.

    if env_path:
        return env_path
    elif os.path.isfile(cwd_path):
        return cwd_path
    else:
        return def_path


def config(path):
    """
    Set CONFIG_PATH, the path to the configuration JSON file,
    to the given path, and return the new path.

    No error or sanity checking is performed.
    """

    global CONFIG_PATH

    CONFIG_PATH = path
    return CONFIG_PATH


def set_config_filename(filename):
    """
    Change the "default" filename for the config file from the
    current value of CONFIG_FILE_NAME to the given filename string.

    The filename is assumed to be a simple filename (no directory
    components).  Partial directory paths might work, or might not.
    No sanity checking is done on the filename.

    This function has no useful effect if CONFIG_PATH is non-None;
    it is used only to compute the CONFIG_PATH, and if the CONFIG_PATH
    has already been set, then CONFIG_PATH will not change as a
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
