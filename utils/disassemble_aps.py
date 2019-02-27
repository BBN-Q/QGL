#!/usr/bin/env python

import numpy as np
import sys
import os.path

from PyQt5.QtWidgets import QMainWindow, QApplication, QWidget, QTableWidget, QTableWidgetItem, QVBoxLayout, QAbstractItemView, QPushButton
from PyQt5.QtGui import QIcon, QColor, QFont
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

colors = {"WFM": QColor(0,200,0),
          "GOTO": QColor(0,100,100),
          "MARKER": QColor(150,150,200),
          "CUSTOM": QColor(200,65,200),
          "WRITEADDR": QColor(245, 105, 65),
          "INVALIDATE": QColor(245, 105, 65),
          "CALL": QColor(65, 205, 245),
          "RET": QColor(65, 205, 245),
          "LOADCMP": QColor(245, 225, 65),
          "MODULATION": QColor(175, 255, 185)}

table_font = QFont("Arial", weight=QFont.Bold)

class MatplotlibWidget(QWidget):
    def __init__(self, I, Q, parent=None):
        super(MatplotlibWidget, self).__init__(parent)
        self.title  = 'Waveform'
        self.left   = 100
        self.top    = 100
        self.width  = 800
        self.height = 600
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        self.figure = Figure()
        self.canvas = FigureCanvasQTAgg(self.figure)

        self.axis = self.figure.add_subplot(111)
        self.axis.plot(I)
        self.axis.plot(Q)
        self.layout = QVBoxLayout(self)
        self.layout.addWidget(self.canvas)
        self.setLayout(self.layout)
        self.canvas.draw()
        self.show()


class DisassemblerApp(QWidget):

    COLUMN_COUNT = 7

    def __init__(self, instructions, waveforms):
        super().__init__()
        self.title  = 'APS2 Disassembled Instructions'
        self.left   = 100
        self.top    = 100
        self.width  = 1000
        self.height = 1200
        self.instructions = instructions
        self.waveforms = waveforms
        self.initUI()
        self.plotters = []

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        self.createTable()
        self.layout = QVBoxLayout()
        self.layout.addWidget(self.tableWidget)
        self.setLayout(self.layout)

        # Show widget
        self.show()

    def createTable(self):
       # Create table
        self.tableWidget = QTableWidget()
        self.tableWidget.setRowCount(len(self.instructions))
        self.tableWidget.setColumnCount(self.COLUMN_COUNT)

        for k, instr in enumerate(self.instructions):
            fields = str(instr).replace(',','').replace(';', '').split(" ")
            if "|" in fields:
                fields.remove("|")
            if fields[0] in colors:
                color = colors[fields[0]]
            else:
                color = None
            for l, f in enumerate(fields):
                text = fields[l]
                if text == "GOTO":
                    btn = QPushButton(self.tableWidget)
                    btn.setText('GOTO')
                    target_row = int(fields[1].split("=")[1])
                    def scroll_to_goto_target(row=target_row, tab=self.tableWidget):
                        tab.scrollToItem(tab.item(row, 0))
                    btn.clicked.connect(scroll_to_goto_target)
                    btn.setStyleSheet("background-color: #006464; color: white; font-weight: bold; text-align: left;")
                    self.tableWidget.setCellWidget(k, l, btn)
                if text == "WFM" and int(fields[4].split("=")[1])==0:
                    # Not a TA pair
                    btn = QPushButton(self.tableWidget)
                    btn.setText(' WFM')
                    addr = int(fields[6].split("=")[1])
                    count = int(fields[5].split("=")[1])
                    def open_plotter(addr=None, I=self.waveforms[0][addr:addr+count], Q=self.waveforms[1][addr:addr+count]):
                        w = MatplotlibWidget(I,Q)
                        self.plotters.append(w)
                    btn.clicked.connect(open_plotter)
                    btn.setStyleSheet("background-color: #00C800; color: white; font-weight: bold; text-align: left;")
                    self.tableWidget.setCellWidget(k, l, btn)
                else:
                    item = QTableWidgetItem(text)
                    item.setFont(table_font)
                    if color:
                        item.setBackground(color)
                    self.tableWidget.setItem(k, l, item)
            if l < self.COLUMN_COUNT-1:
                for j in range(l+1, self.COLUMN_COUNT):
                    item = QTableWidgetItem("")
                    if color:
                        item.setBackground(color)
                    self.tableWidget.setItem(k, j, item)

        self.tableWidget.move(0,0)
        self.tableWidget.setSelectionBehavior(QAbstractItemView.SelectRows)

def get_target_hardware(filename):
    with open(filename, 'rb') as FID:
        target_hw = FID.read(4).decode('utf-8')
    return target_hw

if __name__ == '__main__':
    if len(sys.argv) == 2:
        filename = sys.argv[1]
        try:
            hw = get_target_hardware(filename)
            if hw == "APS2":
                from QGL.drivers.APS2Pattern import read_instructions, read_waveforms
                app = QApplication(sys.argv[:1])
                ex  = DisassemblerApp(read_instructions(sys.argv[1]), read_waveforms(sys.argv[1]))
                sys.exit(app.exec_())
            elif hw == "APS1":
                print("APS1 disassmbler not yet implemented")
            else:
                raise Exception("Unknown hardware target.")
        except:
            print(f"The supplied file {filename} does not appear to be APS1 or APS2 format.")
        
    else:
        print("Expected single file (.aps1 or .aps2 format) as an argument.")

        


