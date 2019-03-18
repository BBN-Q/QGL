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
from ipywidgets import interact, VBox, IntSlider
import ipywidgets as widgets
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


def plot_pulse_files(metafile, time=True, backend='bqplot'):
    '''
    plot_pulse_files(metafile)

    Helper function to plot a list of AWG files. A jupyter slider widget allows choice of sequence number.
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

    line_names, num_seqs, data_dicts = extract_waveforms(fileNames, time=time)
    localname = os.path.split(fileNames[0])[1]
    sequencename = localname.split('-')[0]

    if backend=='matplotlib':
        import matplotlib.pyplot as plt
        def update_plot(seq_ind):
            for line_name in line_names:
                dat = data_dicts[f"{line_name}_{seq_ind}"]
                plt.plot(dat['x'], dat['y'], label=line_name, linewidth=1.0)
        interact(update_plot, seq_ind=IntSlider(min=1,max=num_seqs,step=1,value=1,description="Sequence Number"))

    elif backend=='bqplot':
        from bqplot import DateScale, LinearScale, Axis, Lines, Figure, Tooltip
        from bqplot.colorschemes import CATEGORY10, CATEGORY20
        from ipywidgets import interact, IntSlider
        sx = LinearScale()
        sy = LinearScale(min=-1.0, max=2*len(line_names)-1.0)
        if time:
            ax = Axis(label='Time (ns)', scale=sx)
        else:
            ax = Axis(label="Samples", scale=sx)
        ay = Axis(label='Amplitude', scale=sy, orientation='vertical')

        colors = CATEGORY10 if len(line_names)<10 else CATEGORY20
        lines = []
        tt = Tooltip(fields=['name'], labels=['Channel'])
        x_mult = 1.0e9 if time else 1
        for i, line_name in enumerate(line_names):
            dat = data_dicts[f"{line_name}_1"]
            lines.append(Lines(labels=[line_name], x=x_mult*dat['x'], y=dat['y'], scales={'x': sx, 'y': sy}, 
                                tooltip=tt, animate=False, colors=[colors[i]]))

        slider = IntSlider(min=1, max=num_seqs, step=1, description='Segment', value=1)
        def segment_changed(change):
            for line, line_name in zip(lines, line_names):
                dat = data_dicts[f"{line_name}_{slider.value}"]
                line.x = x_mult*dat['x']
                line.y = dat['y']
        slider.observe(segment_changed, 'value')
        fig = Figure(marks=lines, axes=[ax, ay], title='Waveform Plotter',animation_duration=50)
        return VBox([slider, fig])


def extract_waveforms(fileNames, nameDecorator='', time=False):
    line_names = []
    data_dicts = {}
    num_seqs = 0
    for fileName in sorted(fileNames):

        # Assume a naming convention path/to/file/SequenceName-AWGName.h5
        AWGName = "-".join((os.path.split(os.path.splitext(fileName)[0])[1]).split('-')[1:])
        # Strip any _ suffix
        # if '_' in AWGName:
        #     AWGName = AWGName[:AWGName.index('_')]

        translator = resolve_translator(fileName, translators)
        wfs = translator.read_sequence_file(fileName)
        sample_time = 1.0/translator.SAMPLING_RATE if time else 1

        for (k, seqs) in sorted(wfs.items()):
            if all_zero_seqs(seqs):
                continue
            num_seqs = max(num_seqs, len(seqs))
            line_names.append((AWGName + nameDecorator + '_' + k).replace("-","_"))
            k_ = line_names[-1].replace("-", "_")
            for ct, seq in enumerate(seqs):
                data_dicts[k_ + "_{:d}".format(ct + 1)] = {}

                # Convert from time amplitude pairs to x,y lines with points at start and beginnning to prevent interpolation
                data_dicts[k_ + "_{:d}".format(ct + 1)]["x"] = sample_time * np.tile(
                    np.cumsum([0] + [_[0] for _ in seq]),
                    (2, 1)).flatten(order="F")[1:-1]

                data_dicts[k_ + "_{:d}".format(ct + 1)]["y"] = np.tile(
                    [_[1] for _ in seq],
                    (2, 1)).flatten(order="F") + 2 * (len(line_names) - 1)
    return line_names, num_seqs, data_dicts


# def plot_pulse_files_compare(metafile1, metafile2, time=False):
#     '''
#     plot_pulse_files_compare(fileNames1, fileNames2)

#     Helper function to plot a list of AWG files. A JS slider allows choice of sequence number.
#     '''
#     #If we only go one filename turn it into a list
#     fileNames1 = []
#     fileNames2 = []

#     with open(metafile1, 'r') as FID:
#         meta_info1 = json.load(FID)

#     for el in meta_info1["instruments"].values():
#         # Accomodate seq_file per instrument and per channel
#         if isinstance(el, str):
#             fileNames1.append(el)
#         elif isinstance(el, dict):
#             for file in el.values():
#                 fileNames1.append(file)

#     with open(metafile2, 'r') as FID:
#         meta_info2 = json.load(FID)

#     for el in meta_info2["instruments"].values():
#         # Accomodate seq_file per instrument and per channel
#         if isinstance(el, str):
#             fileNames2.append(el)
#         elif isinstance(el, dict):
#             for file in el.values():
#                 fileNames2.append(file)

#     line_names1, num_seqs1, data_dicts1 = extract_waveforms(fileNames1, 'A', time=time)
#     line_names2, num_seqs2, data_dicts2 = extract_waveforms(fileNames2, 'B', time=time)
#     num_seqs = max(num_seqs1, num_seqs2)
#     data_dicts1.update(data_dicts2)
#     line_names = line_names1 + line_names2

#     localname = os.path.split(fileNames1[0])[1]
#     sequencename = localname.split('-')[0]

#     data = {line_name: ColumnDataSource(data=dat) for line_name,dat in data_dicts1.items()}
#     plot = Figure(title=sequencename, plot_width=1000)
#     plot.background_fill_color = config.plotBackground
#     if config.gridColor:
#         plot.xgrid.grid_line_color = config.gridColor
#         plot.ygrid.grid_line_color = config.gridColor

#     # Colobrewer2 qualitative Set1 (http://colorbrewer2.org)
#     colours = [
#         "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33",
#         "#a65628", "#f781bf", "#999999"
#     ]

#     js_sources = data
#     for ct, k in enumerate(line_names):
#         k_ = k.replace("-", "_")
#         line = plot.line(data_dicts1[k + "_1"]["x"],
#                          data_dicts1[k + "_1"]["y"],
#                          color=colours[ct % len(colours)],
#                          line_width=2,
#                          legend=k.replace("_", "-"))
#         js_sources[k_] = line.data_source

#     code_template = Template("""
#         var seq_num = cb_obj.getv('value');
#         // console.log(seq_num);
#         var all_data = {
#             {% for trace in trace_names %}
#                 '{{trace}}': {{trace}}.getv('data'),
#             {% endfor %}
#         };
#         {% for line in line_names %}
#         {{line}}['data'] = {'x': (all_data['{{line}}_'.concat(seq_num.toString())])['x'],
#                               'y': (all_data['{{line}}_'.concat(seq_num.toString())])['y']};
#         {% endfor %}
#     """)
#     rendered = code_template.render(
#             line_names=line_names,
#             trace_names=list(data.keys())[:-len(line_names)]
#             )

#     callback = CustomJS(
#         args=dict(**js_sources),
#         code=rendered)

#     if num_seqs > 1:
#         slider = Slider(start=1,
#                         end=num_seqs,
#                         value=1,
#                         step=1,
#                         title="Sequence",
#                         callback=callback)

#         layout = column(slider, plot)
#         show(layout)
#     else:
#         show(plot)
