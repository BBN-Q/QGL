
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

import time
import os.path
from importlib import import_module
import pkgutil
from collections import OrderedDict

import numpy as np

from bokeh.layouts import column
from bokeh.models import Slider
from bokeh.plotting import Figure
from bokeh.client import push_session
from bokeh.io import curdoc, output_notebook, output_server

from . import config
from . import ChannelLibrary

from . import drivers

from .Plotting import BokehServerThread, in_notebook

def all_zero_seqs(seqs):
    return all([np.allclose([_[1] for _ in seq], 0) for seq in seqs])

def build_awg_translator_map():
    translators_map = {}
    translators = [_[1] for _ in pkgutil.walk_packages(drivers.__path__)]
    for translator in translators:
        module = import_module('QGL.drivers.' + translator)
        ext = module.get_seq_file_extension()
        if ext in translators_map:
            translators_map[ext].append(module)
        else:
            translators_map[ext] = [module]
    return translators_map

# static translator map
translators = build_awg_translator_map()

def resolve_translator(filename, translators):
    ext = os.path.splitext(filename)[1]
    if ext not in translators:
        raise NameError("No translator found to open the given file %s", filename)
    if len(translators[ext]) == 1:
        return translators[ext][0]
    for t in translators[ext]:
        if t.is_compatible_file(filename):
            return t
    raise NameError("No translator found to open the given file %s", filename)

def plot_pulse_files(h5_files):
    """
    Plots a list of h5 sequence files with an interactive slider.
    """

    bokeh_thread = BokehServerThread()
    bokeh_thread.start()
    time.sleep(1) #make sure server is finish launching
    output_server()
    curdoc().clear()
    # output_notebook()

    #If we only got one filename turn it into a list
    if isinstance(h5_files, str):
        h5_files = [h5_files]

    wfs = extract_waveforms(h5_files)
    num_seqs = max([len(_) for _ in wfs.values()])

    filename = os.path.split(h5_files[0])[1]
    seq_name = filename.split('-')[0]

    fig = Figure(title=seq_name, plot_width=800)
    fig.background_fill_color = config.plotBackground
    if config.gridColor:
        fig.xgrid.grid_line_color = config.gridColor
        fig.ygrid.grid_line_color = config.gridColor

    num_lines = len(wfs.keys())
    #for some reason the qualitative maps aren't in bokeh.palettes
    # see https://github.com/bokeh/bokeh/issues/4758
    brewer_set3 = ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f']
    colours = list(np.resize(brewer_set3, num_lines))

    lines = []
    for colour, (k,v) in zip(colours, wfs.items()):
        lines.append( fig.line(v[0]["x"], v[0]["y"], color=colour, legend=k, line_width=2) )

    def callback(attr, old, new):
        for line, wf in zip(lines, wfs.values()):
            new_xy = {"x":wf[new-1]["x"], "y":wf[new-1]["y"]}
            line.data_source.data = new_xy

    slider = Slider(start=1, end=num_seqs, value=1, step=1, title="Sequence")
    slider.on_change("value", callback)

    session = push_session(curdoc())
    session.show(column([slider, fig]))
    session.loop_until_closed()

def extract_waveforms(h5_files, nameDecorator=''):
    """
    Extracts a dictionary mapping channel names to lists of dictionaries (keys "x" and "y") lines

    Lines are shifted along the y axis to prevent overlap
    """
    wfs = OrderedDict()

    y_shift = 0
    for file_name in sorted(h5_files):
        # Assume a naming convention path/to/file/SequenceName-AWGName.h5
        AWGName = (os.path.split(os.path.splitext(file_name)[0])[1]).split('-')[1]
        # Strip any _ suffix
        if '_' in AWGName:
            AWGName = AWGName[:AWGName.index('_')]

        translator = resolve_translator(file_name, translators)
        ta_pairs = translator.read_sequence_file(file_name)

        for (k,seqs) in sorted(ta_pairs.items()):
            if all_zero_seqs(seqs):
                continue

            line_name = AWGName + nameDecorator + '-' + k
            y_shift += 2
            wfs[line_name] = []
            for ct,seq in enumerate(seqs):
                wfs[line_name].append(
                    {
                        # Convert from time amplitude pairs to x,y lines with points at start and beginnning to prevent interpolation
                        "x": np.tile( np.cumsum([0] + [_[0] for _ in seq]), (2,1)).flatten(order="F")[1:-1],
                        "y": np.tile( [_[1] for _ in seq], (2,1)).flatten(order="F") + y_shift
                    }
                )
    return wfs

def plot_pulse_files_compare(fileNames1, fileNames2):
    '''
    plot_pulse_files_compare(fileNames1, fileNames2)

    Helper function to plot a list of AWG files. A JS slider allows choice of sequence number.
    '''
    #If we only go one filename turn it into a list
    if isinstance(fileNames1, str):
        fileNames1 = [fileNames1]
    if isinstance(fileNames2, str):
        fileNames2 = [fileNames2]

    dataDict = {}

    lineNames1, num_seqs1 = extract_waveforms(dataDict, fileNames1, 'A')
    lineNames2, num_seqs2 = extract_waveforms(dataDict, fileNames2, 'B')
    num_seqs = max(num_seqs1, num_seqs2)

    localname = os.path.split(fileNames1[0])[1]
    seq_name = localname.split('-')[0]

    all_data = ColumnDataSource(data=dataDict)
    plot = Figure(title=seq_name, plot_width=1000)
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
    for ct,k in enumerate(lineNames1 + lineNames2):
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

    callback = CustomJS(args=js_sources, code=code_template.render(lineNames=[l.replace("-","_") for l in lineNames1+lineNames2]))

    slider = Slider(start=1, end=num_seqs, value=1, step=1, title="Sequence", callback=callback)

    layout = vform(slider, plot)

    show(layout)
