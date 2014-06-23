# -*- coding: utf-8 -*-

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import matplotlib.pyplot as plt

import random

class CentralWidget(QWidget):

    def __init__(self):
        super(CentralWidget, self).__init__()
 
        # Widgets
        self.lineEditDirectory = QLineEdit(self)
        self.pushButtonSelect = QPushButton("...", self)
        self.comboBoxParamTag = QComboBox(self)
        self.comboBoxDataTag = QComboBox(self)
        self.listParametersX = QListWidget(self)
        self.listParametersY = QListWidget(self)
        self.pushButtonPlot = QPushButton("Make plot", self)

        # Connections





        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        self.toolbar = NavigationToolbar(self.canvas, self)

        # set the layout

        layoutL = QGridLayout()
        layoutL.setSpacing(10)
        layoutL.addWidget(self.lineEditDirectory, 0, 0, 1, 2)
        layoutL.addWidget(self.pushButtonSelect, 0, 2)
        layoutL.addWidget(self.comboBoxParamTag, 1, 0)
        layoutL.addWidget(self.comboBoxDataTag, 1, 1)
        layoutL.addWidget(self.listParametersX, 2, 0)
        layoutL.addWidget(self.listParametersY, 2, 1)
        layoutL.addWidget(self.pushButtonPlot, 3, 2)
        wl = QWidget(self)
        wl.setLayout(layoutL)

        layoutR = QVBoxLayout()
        layoutR.addWidget(self.toolbar)
        layoutR.addWidget(self.canvas)
        wr = QWidget(self)
        wr.setLayout(layoutR)

        splitter = QSplitter(Qt.Horizontal, self)
        splitter.addWidget(wl)
        splitter.addWidget(wr)
        
        layout = QVBoxLayout()
        layout.addWidget(splitter)
        self.setLayout(layout)

        # test with random data
        self.test()



    def select(self):
        pass

    def chooseParamTag(self, strParamTag):
        pass
    
    def chooseDataTag(self, strDataTag):
        pass
    
    def doubleClickedParamX(self, item):
        pass

    def doubleClickedParamY(self, item):
        pass


    def test(self):

        # random data
        data = [random.random() for i in range(10)]

        # create an axis
        ax = self.figure.add_subplot(111)

        # discards the old graph
        ax.hold(False)

        # plot data
        ax.plot(data, '*-')

        # refresh canvas
        self.canvas.draw()

