'''
Created on Feb 18, 2015

@author: cryan@bbn.com and bjohnson@bbn.com

Common plotting tools.

Copyright 2012-2015 Raytheon BBN Technologies

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

import os.path, uuid, tempfile
import bokeh.plotting as bk

#Effective global of whether we are interactive plottinn in a notebook or to a
#static html file
_in_ipynb = False

def output_notebook():
    global _in_ipynb
    bk.output_notebook()
    _in_ipynb = True

def output_file():
    global _in_ipynb
    bk.output_file(os.path.join(tempfile.gettempdir(), str(uuid.uuid4()) + ".html"))
    _in_ipynb = False

def in_ipynb():
    return _in_ipynb

#Default to output_file
output_file()
