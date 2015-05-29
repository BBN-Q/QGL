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

def build_waveforms(seq):
    from Compiler import compile_sequence
    #compile
    linkList, wfLib = compile_sequence(seq)

    # build a concatenated waveform for each channel
    channels = linkList.keys()
    concatShapes = {q: np.array([0], dtype=np.complex128) for q in channels}
    for q in channels:
        # TODO: deal with repeated sections
        for entry in filter(lambda x: isinstance(x, Compiler.LLWaveform), linkList[q]):
            if entry.isTimeAmp:
                concatShapes[q] = np.append(concatShapes[q], wfLib[q][entry.key][0]*np.ones(entry.length))
            else:
                concatShapes[q] = np.append(concatShapes[q], wfLib[q][entry.key])
    # add an extra zero to make things look more normal
    for q in channels:
        concatShapes[q] = np.append(concatShapes[q], 0)  
    return concatShapes


def plot_waveforms(waveforms, figTitle = ''):
    AWGFreq = 1.2e9 # TODO: make this independent of sampling rate

    channels = waveforms.keys()
     # plot
    plots = []
    for (ct,chan) in enumerate(channels):
        fig = bk.figure(title=figTitle + repr(chan), plot_width=800, plot_height=350, y_range=[-1.05, 1.05])
        waveformToPlot = waveforms[chan]
        xpts = np.linspace(0,len(waveformToPlot)/AWGFreq/1e-6,len(waveformToPlot))
        fig.line(xpts, np.real(waveformToPlot), color='red')
        fig.line(xpts, np.imag(waveformToPlot), color='blue')
        plots.append(fig)
    bk.show(bk.VBox(*plots))

def show(seq):
    from Compiler import build_waveforms
    waveforms = build_waveforms(seq)
    plot_waveforms(waveforms)
