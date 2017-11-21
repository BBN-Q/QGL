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

import QGL
import os
import unittest

class ConfigTest(unittest.TestCase):
    """
    Tests the config.py code.
    """

    def setUp(self):
        pass

    def test_env(self):
        """
        Test that if the BBN_MEAS_FILE environment variable is set, it is used
        """

        # Write a bare-bones yml config file
        meas_name = '/tmp/qgl_test_1.yml'
        meas_txt  = 'config:\n'
        meas_txt += '    AWGDir: /foo/bar/xyz\n'
        fout = open(meas_name, 'w')
        fout.write(meas_txt)
        fout.close()

        orig_env = os.getenv('BBN_MEAS_FILE')
        os.environ['BBN_MEAS_FILE'] = meas_name
        QGL.config.load_config() 
        os.environ['BBN_MEAS_FILE'] = orig_env
        assert QGL.config.meas_file == meas_name
        assert QGL.config.AWGDir == "/foo/bar/xyz"

    def test_override1(self):
        """
        Tests manually supplying a different config file when instantiating the channel library.
        """

        # Write a bare-bones yml config file
        meas_name = '/tmp/qgl_test_2.yml'
        meas_txt  = 'config:\n'
        meas_txt += '    AWGDir: /foo/bar/abc\n'
        fout = open(meas_name, 'w')
        fout.write(meas_txt)
        fout.close()

        cl = QGL.config.load_config(meas_name) 
        assert QGL.config.meas_file == meas_name
        assert QGL.config.AWGDir == "/foo/bar/abc"

if __name__ == "__main__":
    unittest.main()
