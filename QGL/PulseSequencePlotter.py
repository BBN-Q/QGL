
'''
Created on Jan 8, 2012

Large refactor to notebook Bokeh plotting 5 Jan, 2015

@author: cryan@bbn.com

A simple interactive tool plotting pulse sequences for visual inspection

Copyright 2012,2013,2014,2015 Raytheon BBN Technologies

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

import os.path
from importlib import import_module
from bokeh.io import vform
from bokeh.models import CustomJS, ColumnDataSource, Slider
from bokeh.plotting import Figure, show

from jinja2 import Template
import numpy as np

from . import config
from . import ChannelLibrary

def all_zero_seqs(seqs):
    return all([np.all(seq[1] == 0) for seq in seqs])

def build_awg_translator_map():
    translators = {}
    for chan in ChannelLibrary.channelLib.values():
        if (hasattr(chan, 'AWG') and chan.AWG and
            hasattr(chan, 'translator') and chan.translator):
            translators[chan.AWG] = import_module('QGL.drivers.'+chan.translator)
    return translators

def plot_pulse_files(fileNames):
    '''
    plot_pulse_files(fileNames)

    Helper function to plot a list of AWG files. A JS slider allows choice of sequence number.
    '''
    #If we only go one filename turn it into a list
    if isinstance(fileNames, str):
        fileNames = [fileNames]

    wfs = {}
    dataDict = {}
    lineNames = []
    title = ""
    translators = build_awg_translator_map()
    num_seqs = 0
    for fileName in sorted(fileNames):

        #Assume a naming convention path/to/file/SequenceName-AWGName.h5
        AWGName = (os.path.split(os.path.splitext(fileName)[0])[1]).split('-')[1]
        #Strip any _ suffix
        if '_' in AWGName:
            AWGName = AWGName[:AWGName.index('_')]

        title += os.path.split(os.path.splitext(fileName)[0])[1] + "; "

        wfs[AWGName] = translators[AWGName].read_sequence_file(fileName)

        for (k,seqs) in sorted(wfs[AWGName].items()):
            if not all_zero_seqs(seqs):
                num_seqs = max(num_seqs, len(seqs))
                lineNames.append(AWGName + '-' + k)
                k_ = lineNames[-1].replace("-", "_")
                for ct,seq in enumerate(seqs):
                    #Add in end points for each timestep to give hold rather than interpolation effect
                    dataDict[k_+"_x_{:d}".format(ct+1)] = np.tile(seqs[ct][0], (2,1)).flatten(order="F")[1:]
                    dataDict[k_+"_y_{:d}".format(ct+1)] = np.tile(seqs[ct][1], (2,1)).flatten(order="F")[0:-1] + 2*(len(lineNames)-1)

    #Remove trailing semicolon from title
    title = title[:-2]

    all_data = ColumnDataSource(data=dataDict)
    plot = Figure(title=title, plot_width=1000)
    plot.background_fill_color = config.plotBackground
    if config.gridColor:
        plot.xgrid.grid_line_color = config.gridColor
        plot.ygrid.grid_line_color = config.gridColor

    # Colobrewer2 qualitative Set1 (http://colorbrewer2.org)
    colours = [
        "#e41a1c",
        "#377eb8",
        "#4daf4a",
        "#984ea3",
        "#ff7f00",
        "#ffff33",
        "#a65628",
        "#f781bf",
        "#999999"
    ]

    js_sources = {}
    js_sources["all_data"] = all_data
    for ct,k in enumerate(lineNames):
        k_ = k.replace("-", "_")
        line = plot.line(dataDict[k_+"_x_1"], dataDict[k_+"_y_1"], color=colours[ct%len(colours)], line_width=2, legend=k)
        js_sources[k_] = line.data_source

    code_template = Template("""
        var seq_num = cb_obj.get('value');
        console.log(seq_num)
        var all_data = all_data.get('data');
        {% for line in lineNames %}
        {{line}}.set('data', {'x':all_data['{{line}}_x_'.concat(seq_num.toString())], 'y':all_data['{{line}}_y_'.concat(seq_num.toString())]} );
        {% endfor %}
        console.log("Got here!")
    """)

    callback = CustomJS(args=js_sources, code=code_template.render(lineNames=[l.replace("-","_") for l in lineNames]))

    slider = Slider(start=1, end=num_seqs, value=1, step=1, title="Sequence", callback=callback)

    layout = vform(slider, plot)

    show(layout)
