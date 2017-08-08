# -*- coding: UTF-8 -*-
'''
Created on August 2, 2017

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

Tests for the changes to config.py and qgl_config_loc.py
'''

import os
import subprocess
import unittest

class ConfigTest(unittest.TestCase):
    """
    Tests the config.py and qgl_config_loc.py, focusing on the latter.

    Because the relevant code is executed at import time, the different
    effects are awkward/impossible to test within a single program.
    Therefore most of the test cases are short programs that are run
    in a separate Python process, and the output of these programs is
    checked against the expected output.

    In order to keep this test self-contained, the programs are
    embedded in the tests, rather than represented in separate files.
    """

    def setUp(self):
        pass

    def test_default(self):
        """
        Test that if config.py is imported in the 'default' manner,
        the default behavior is observed
        """

        qgl_dir = os.path.dirname(os.path.dirname(__file__))

        os.putenv('PYTHONPATH', qgl_dir)

        progname = '/tmp/config_loc_test.py'

        progtext = 'import qgl_config_loc\n'
        progtext += 'import QGL.config\n'

        fout = open(progname, 'w')
        fout.write(progtext)
        fout.close()

        proc = subprocess.Popen(
                ['python', progname],
                universal_newlines=True, stdout=subprocess.PIPE)
        out_data, err_data = proc.communicate()

        lines = out_data.split('\n')
        expected = 'Note: using QGLCfgFile [%s/QGL/config.json]' % qgl_dir
        assert lines[0] == expected

    def test_env(self):
        """
        Test that if the QGLCFGFILE environment variable is set, it is used
        """

        qgl_dir = os.path.dirname(os.path.dirname(__file__))

        os.putenv('PYTHONPATH', qgl_dir)
        os.putenv('QGLCFGFILE', '/foo/bar/qux')

        progname = '/tmp/config_loc_test.py'

        progtext = 'import qgl_config_loc\n'
        progtext += 'import QGL.config\n'

        fout = open(progname, 'w')
        fout.write(progtext)
        fout.close()

        try:
            proc = subprocess.Popen(
                    ['python', progname],
                    universal_newlines=True,
                    stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            out_data, err_data = proc.communicate()
        except BaseException as exc:
            assert False

        os.unsetenv('QGLCFGFILE')

        lines = out_data.split('\n')
        assert lines[0] == 'Note: using QGLCfgFile [/foo/bar/qux]'

    def test_override1(self):
        """
        Test that if the config path is overridden, the new value is used

        The first test just uses a default environment.
        """

        qgl_dir = os.path.dirname(os.path.dirname(__file__))

        os.putenv('PYTHONPATH', qgl_dir)

        progname = '/tmp/config_loc_test.py'

        progtext = 'import qgl_config_loc\n'
        progtext += 'qgl_config_loc.config(\'/a/b/c/d\')\n'
        progtext += 'import QGL.config\n'

        fout = open(progname, 'w')
        fout.write(progtext)
        fout.close()

        try:
            proc = subprocess.Popen(
                    ['python', progname],
                    universal_newlines=True,
                    stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            out_data, err_data = proc.communicate()
        except BaseException as exc:
            assert False

        lines = out_data.split('\n')
        assert lines[0] == 'Note: using QGLCfgFile [/a/b/c/d]'


if __name__ == "__main__":
    unittest.main()
