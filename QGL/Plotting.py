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
import numpy as np

import threading
import subprocess
import psutil
import os
import sys

from . import config


def output_notebook():
    bk.output_notebook()


def output_file():
    bk.output_file(os.path.join(tempfile.gettempdir(), str(uuid.uuid4()) +
                                ".html"))

# Default to output_file
# Commented out because recent versions of Bokeh remember the state of
# output_notebook/output_file on a per module basis.

# output_file()


def build_waveforms(seq):
    # import here to avoid circular imports
    from . import Compiler, PulseSequencer, BlockLabel, ControlFlow
    wires = Compiler.compile_sequence(seq)

    # build a concatenated waveform for each channel
    channels = wires.keys()
    concatShapes = {q: np.array([0], dtype=np.complex128) for q in channels}
    for q in channels:
        frame = 0
        repeat = 0
        label_idx = {}
        for entry in wires[q]:
            if isinstance(entry, PulseSequencer.Pulse):
                shape = np.exp(1j *
                               (frame + entry.phase)) * entry.amp * entry.shape
                frame += entry.frameChange
                concatShapes[q] = np.append(concatShapes[q], shape)
            elif isinstance(entry, BlockLabel.BlockLabel):
                label_idx[entry] = len(concatShapes[q]) + 1
            elif isinstance(entry, ControlFlow.LoadRepeat):
                repeat = entry.value - 1
            elif isinstance(entry, ControlFlow.Repeat):
                concatShapes[q] = np.append(concatShapes[q], np.tile(
                    concatShapes[q][label_idx[entry.target]:], repeat))

    # add an extra zero to make things look more normal
    for q in channels:
        concatShapes[q] = np.append(concatShapes[q], 0)
    return concatShapes


def plot_waveforms(waveforms, figTitle=''):
    channels = waveforms.keys()
    # plot
    plots = []
    for (ct, chan) in enumerate(channels):
        fig = bk.figure(title=figTitle + repr(chan),
                        plot_width=800,
                        plot_height=350,
                        y_range=[-1.05, 1.05])
        fig.background_fill_color = config.plotBackground
        if config.gridColor:
            fig.xgrid.grid_line_color = config.gridColor
            fig.ygrid.grid_line_color = config.gridColor
        waveformToPlot = waveforms[chan]
        xpts = np.linspace(0, len(waveformToPlot) / chan.physChan.samplingRate
                           / 1e-6, len(waveformToPlot))
        fig.line(xpts, np.real(waveformToPlot), color='red')
        fig.line(xpts, np.imag(waveformToPlot), color='blue')
        plots.append(fig)
    bk.show(bk.vplot(*plots))


def show(seq):
    waveforms = build_waveforms(seq)
    plot_waveforms(waveforms)

def in_notebook():
    return True
    # From http://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook
    try:
        cfg = get_ipython().config
        if cfg['IPKernelApp']['parent_appname'] == 'ipython-notebook':
            return True
        else:
            return False
    except NameError:
        return False

class BokehServerThread(threading.Thread):
    def __init__(self):
        super(BokehServerThread, self).__init__()
        self.daemon = True
        self.run_in_notebook = in_notebook()

    def __del__(self):
        self.join()

    def __del__(self):
        self.join()

    def run(self):
        args = ["bokeh", "serve"]
        if self.run_in_notebook:
            args.append("--allow-websocket-origin=localhost:8888")
        self.p = subprocess.Popen(args, env=os.environ.copy())

    def join(self, timeout=None):
        if self.p:
            print("Killing bokeh server thread {}".format(self.p.pid))
            for child_proc in psutil.Process(self.p.pid).children():
                print("Killing child process {}".format(child_proc.pid))
                child_proc.kill()
            self.p.kill()
            self.p = None
            super(BokehServerThread, self).join(timeout=timeout)
