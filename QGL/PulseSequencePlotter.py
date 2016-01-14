
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
import bokeh.plotting as bk
from bokeh.embed import notebook_div
from Plotting import in_ipynb

from IPython.html import widgets
from IPython.display import display

import numpy as np

import uuid, tempfile

import config

def all_zero_seqs(seqs):
    return all([np.all(seq == 0) for seq in seqs])

def plot_pulse_files(fileNames, firstSeqNum=0):
    '''
    plot_pulse_files(fileNames, firstSeqNum=0)

    Helper function to plot a list of AWG files.  In an iPython notebook the plots will be in line with
    dynamic updating.  For iPython consoles a static html file will be generated with the firstSeqNum.
    '''
    #If we only go one filename turn it into a list
    if isinstance(fileNames, str):
        fileNames = [fileNames]

    wfs = {}
    dataDict = {}
    lineNames = []
    title = ""

    import Libraries # prevent circular import

    for fileName in sorted(fileNames):

        #Assume a naming convention path/to/file/SequenceName-AWGName.h5
        AWGName = (os.path.split(os.path.splitext(fileName)[0])[1]).split('-')[1]
        #Strip any _ suffix
        if '_' in AWGName:
            AWGName = AWGName[:AWGName.index('_')]

        title += os.path.split(os.path.splitext(fileName)[0])[1] + "; "

        wfs[AWGName] = Libraries.instrumentLib[AWGName].read_sequence_file(fileName)

        for (k,seqs) in sorted(wfs[AWGName].items()):
            if not all_zero_seqs(seqs):
                lineNames.append(AWGName + '-' + k)
                dataDict[lineNames[-1] + "_x"] = np.arange(len(seqs[firstSeqNum]))
                dataDict[lineNames[-1]] = seqs[firstSeqNum] + 2*(len(lineNames)-1)

    #Remove trailing semicolon from title
    title = title[:-2]

    source = bk.ColumnDataSource(data=dataDict)
    figH = bk.figure(title=title, plot_width=1000, y_range=(-1,len(dataDict)+1))
    figH.background_fill = config.plotBackground
    if config.gridColor:
        figH.xgrid.grid_line_color = config.gridColor
        figH.ygrid.grid_line_color = config.gridColor

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

    for ct,k in enumerate(lineNames):
        figH.line(k+"_x", k, source=source, color=colours[ct%len(colours)], line_width=2, legend=k)

    if in_ipynb():
        #Setup inline plotting with slider updating of data

        def update_plot(_, seqNum):
            for ct,k in enumerate(lineNames):
                AWGName, chName = k.split('-')
                source.data[k+"_x"] = np.arange(len(wfs[AWGName][chName][seqNum-1]))
                source.data[k] = wfs[AWGName][chName][seqNum-1] + 2*ct
            source.push_notebook()

        #widgets.interact(update_plot, seqNum=(1, len(wfs[AWGName]["ch1"])), div=widgets.HTMLWidget(value=notebook_div(figH)))
        if 'ch1' in wfs[AWGName].keys():
            chkey = "ch1"
        else:
            chkey = "ch3"
        slider = widgets.IntSlider(value=firstSeqNum+1, min=1, max=len(wfs[AWGName][chkey]), step=1, description="Sequence (of {}):".format(len(seqs)))
        slider.on_trait_change(update_plot, 'value')
        plotBox = widgets.HTML(value=notebook_div(figH))
        appBox = widgets.Box()
        appBox.children = [slider, plotBox]
        display(appBox)

    else:
        #Otherwise dump to a static file
        bk.show(figH)
