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
import json
from importlib import import_module
from bokeh.layouts import column
from bokeh.models import CustomJS, ColumnDataSource, Slider
from bokeh.plotting import Figure, show

from jinja2 import Template
import numpy as np

from . import config
from . import drivers
import pkgutil


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
        raise NameError("No translator found to open the given file %s",
                        filename)
    if len(translators[ext]) == 1:
        return translators[ext][0]
    for t in translators[ext]:
        if t.is_compatible_file(filename):
            return t
    raise NameError("No translator found to open the given file %s", filename)


def plot_pulse_files(metafile, time=False):
    '''
    plot_pulse_files(metafile)

    Helper function to plot a list of AWG files. A JS slider allows choice of sequence number.
    '''
    #If we only go one filename turn it into a list
    
    with open(metafile, 'r') as FID:
        meta_info = json.load(FID)
    fileNames = []
    for el in meta_info["instruments"].values():
        # Accomodate seq_file per instrument and per channel
        if isinstance(el, str):
            fileNames.append(el)
        elif isinstance(el, dict):
            for file in el.values():
                fileNames.append(file)

    dataDict = {}
    lineNames, num_seqs = extract_waveforms(dataDict, fileNames, time=time)

    localname = os.path.split(fileNames[0])[1]
    sequencename = localname.split('-')[0]

    all_data = ColumnDataSource(data=dataDict)
    plot = Figure(title=sequencename, plot_width=1000)
    plot.background_fill_color = config.plotBackground
    if config.gridColor:
        plot.xgrid.grid_line_color = config.gridColor
        plot.ygrid.grid_line_color = config.gridColor

    # Colobrewer2 qualitative Set1 (http://colorbrewer2.org)
    colours = [
        "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33",
        "#a65628", "#f781bf", "#999999"
    ]

    js_sources = {}
    js_sources["all_data"] = all_data
    for ct, k in enumerate(lineNames):
        k_ = k.replace("-", "_")
        line = plot.line(dataDict[k_ + "_x_1"],
                         dataDict[k_ + "_y_1"],
                         color=colours[ct % len(colours)],
                         line_width=2,
                         legend=k)
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

    callback = CustomJS(
        args=js_sources,
        code=code_template.render(
            lineNames=[l.replace("-", "_") for l in lineNames]))

    slider = Slider(start=1,
                    end=num_seqs,
                    value=1,
                    step=1,
                    title="Sequence",
                    callback=callback)

    layout = column(slider, plot)

    show(layout)


def extract_waveforms(dataDict, fileNames, nameDecorator='', time=False):
    lineNames = []
    num_seqs = 0
    for fileName in sorted(fileNames):

        # Assume a naming convention path/to/file/SequenceName-AWGName.h5
        AWGName = "-".join((os.path.split(os.path.splitext(fileName)[0])[1]).split('-')[1:])
        # Strip any _ suffix
        if '_' in AWGName:
            AWGName = AWGName[:AWGName.index('_')]

        translator = resolve_translator(fileName, translators)
        wfs = translator.read_sequence_file(fileName)
        sample_time = 1.0/translator.SAMPLING_RATE if time else 1

        for (k, seqs) in sorted(wfs.items()):
            if all_zero_seqs(seqs):
                continue
            num_seqs = max(num_seqs, len(seqs))
            lineNames.append(AWGName + nameDecorator + '-' + k)
            k_ = lineNames[-1].replace("-", "_")
            for ct, seq in enumerate(seqs):
                # Convert from time amplitude pairs to x,y lines with points at start and beginnning to prevent interpolation
                dataDict[k_ + "_x_{:d}".format(ct + 1)] = sample_time * np.tile(
                    np.cumsum([0] + [_[0] for _ in seq]),
                    (2, 1)).flatten(order="F")[1:-1]
                dataDict[k_ + "_y_{:d}".format(ct + 1)] = np.tile(
                    [_[1] for _ in seq],
                    (2, 1)).flatten(order="F") + 2 * (len(lineNames) - 1)
    return lineNames, num_seqs


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
    sequencename = localname.split('-')[0]

    all_data = ColumnDataSource(data=dataDict)
    plot = Figure(title=sequencename, plot_width=1000)
    plot.background_fill_color = config.plotBackground
    if config.gridColor:
        plot.xgrid.grid_line_color = config.gridColor
        plot.ygrid.grid_line_color = config.gridColor

    # Colobrewer2 qualitative Set1 (http://colorbrewer2.org)
    colours = [
        "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33",
        "#a65628", "#f781bf", "#999999"
    ]

    js_sources = {}
    js_sources["all_data"] = all_data
    for ct, k in enumerate(lineNames1 + lineNames2):
        k_ = k.replace("-", "_")
        line = plot.line(dataDict[k_ + "_x_1"],
                         dataDict[k_ + "_y_1"],
                         color=colours[ct % len(colours)],
                         line_width=2,
                         legend=k)
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

    callback = CustomJS(
        args=js_sources,
        code=code_template.render(
            lineNames=[l.replace("-", "_") for l in lineNames1 + lineNames2]))

    slider = Slider(start=1,
                    end=num_seqs,
                    value=1,
                    step=1,
                    title="Sequence",
                    callback=callback)

    layout = column(slider, plot)

    show(layout)
