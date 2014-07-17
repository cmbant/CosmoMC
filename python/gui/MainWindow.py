# -*- coding: utf-8 -*-

import os
import glob
import signal
import random

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import matplotlib.pyplot as plt

import GetDistPlots

class MainWindow(QMainWindow):

    def __init__(self):

        super(MainWindow, self).__init__()

        # 
        self.setAttribute(Qt.WA_DeleteOnClose)
 
        self.createActions()
        self.createMenus()
        #self.createToolBars()
        self.createStatusBar()
        self.createDockWidgets()
        self.createCentralWidget()

        self.setWindowTitle("GetDist GUI")
        self.resize(800, 600)

        # Allow to shutdown the GUI with Ctrl+C
        signal.signal(signal.SIGINT, signal.SIG_DFL)

        # 
        self.rootdir = ""


    def createActions(self):

        self.exportAct = QAction(
            "&Export", self,
            statusTip="Export image as PDF", triggered=self.export)

        self.scriptAct = QAction(
            "Script", self,
            statusTip="Export commands to script", triggered=self.script)

        self.exitAct = QAction(QIcon("./images/application_exit.png"), 
            "E&xit", self, shortcut="Ctrl+Q",
            statusTip="Exit application",
            triggered=qApp.closeAllWindows)
 
        self.aboutAct = QAction(
            "&About", self,
            statusTip="Show About box",
            triggered=self.about)

    def createMenus(self):

        self.fileMenu = self.menuBar().addMenu("&File")
        self.fileMenu.addAction(self.exportAct)
        self.separatorAct = self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.scriptAct)
        self.separatorAct = self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.exitAct)
 
        self.menuBar().addSeparator()
 
        self.helpMenu = self.menuBar().addMenu("&Help")
        self.helpMenu.addAction(self.aboutAct)

    def createToolBars(self):
        self.fileToolBar = self.addToolBar("File")
        self.fileToolBar.addAction(self.exportAct)
        self.fileToolBar.addAction(self.exitAct)
 
    def createStatusBar(self):
        self.statusBar().showMessage("Ready", 2000)

    def createDockWidgets(self):

        # Top: select Widget
        dockTop = QDockWidget(self.tr("Filters"), self)
        dockTop.setAllowedAreas(Qt.LeftDockWidgetArea)

        self.selectWidget = QWidget(dockTop)

        self.lineEditDirectory = QLineEdit(self.selectWidget)
        self.lineEditDirectory.clear()
        self.lineEditDirectory.setReadOnly(True)

        self.pushButtonSelect = QPushButton("...", self.selectWidget)
        self.pushButtonSelect.setToolTip("Choose root directory")
        self.connect(self.pushButtonSelect, SIGNAL("clicked()"), self.selectRootDir)

        self.comboBoxParamTag = QComboBox(self.selectWidget)
        self.comboBoxParamTag.clear()
        self.connect(self.comboBoxParamTag, SIGNAL("activated(const QString&)"), self.setParamTag)

        self.comboBoxDataTag = QComboBox(self.selectWidget)
        self.comboBoxDataTag.clear()
        self.connect(self.comboBoxDataTag, SIGNAL("activated(const QString&)"), self.setDataTag)

        self.listParametersX = QListWidget(self.selectWidget)
        self.listParametersX.clear()
        self.connect(self.listParametersX, SIGNAL("itemClicked()"), self.clickParamX)

        self.listParametersY = QListWidget(self.selectWidget)
        self.listParametersY.clear()
        self.connect(self.listParametersY, SIGNAL("itemClicked()"), self.clickParamY)

        self.selectAllX = QCheckBox("Select All", self.selectWidget)
        self.selectAllX.setCheckState(Qt.Unchecked)
        self.connect(self.selectAllX, SIGNAL("clicked()"), self.statusSelectAllX)

        self.selectAllY = QCheckBox("Select All", self.selectWidget)
        self.selectAllY.setCheckState(Qt.Unchecked)
        self.connect(self.selectAllY, SIGNAL("clicked()"), self.statusSelectAllY)

        self.pushButtonPlot = QPushButton("Make plot", self.selectWidget)
        self.connect(self.pushButtonPlot, SIGNAL("clicked()"), self.plotData)

        # Graphic Layout
        layoutTop = QGridLayout()
        layoutTop.setSpacing(5)
        layoutTop.addWidget(self.lineEditDirectory, 0, 0, 1, 3)
        layoutTop.addWidget(self.pushButtonSelect,  0, 3, 1, 1)
        layoutTop.addWidget(self.comboBoxParamTag,  1, 0, 1, 2)
        layoutTop.addWidget(self.comboBoxDataTag,   1, 2, 1, 2)
        layoutTop.addWidget(self.listParametersX,   2, 0, 1, 2)
        layoutTop.addWidget(self.listParametersY,   2, 2, 1, 2)
        layoutTop.addWidget(self.selectAllX,        3, 0, 1, 2)
        layoutTop.addWidget(self.selectAllY,        3, 2, 1, 2)
        layoutTop.addWidget(self.pushButtonPlot,    4, 0, 1, 4)
        self.selectWidget.setLayout(layoutTop)
        
        # Minimum size for initial resize of dock widget
        #self.selectWidget.setMinimumSize(250, 250)
        self.selectWidget.setMaximumSize(300, 300)

        dockTop.setWidget(self.selectWidget)
        self.addDockWidget(Qt.LeftDockWidgetArea, dockTop)


        # Bottom: Tab Widget  
        dockBot = QDockWidget(self.tr("Data Files"), self)
        dockBot.setAllowedAreas(Qt.LeftDockWidgetArea)

        self.dataWidget = QWidget(dockBot)

        self.newButton = QPushButton(QIcon("images/file_new.png"), "New", self.dataWidget)
        self.connect(self.newButton, SIGNAL("clicked()"), self.openNewFile)

        spacer = QSpacerItem(10, 10, QSizePolicy.Expanding, QSizePolicy.Ignored)

        self.tabWidget = QTabWidget(dockBot)
        self.tabWidget.setMovable(False)
        self.tabWidget.setTabsClosable(True)
        self.connect(self.tabWidget, SIGNAL("tabCloseRequested(int)"), self.closeTab)

        # Graphic Layout
        layoutBot = QGridLayout()
        layoutBot.setSpacing(5)
        layoutBot.addWidget(self.newButton, 0, 0, 1, 1)
        layoutBot.addItem(spacer, 0, 1, 1, 1)
        layoutBot.addWidget(self.tabWidget, 1, 0, 1, 2)
        self.dataWidget.setLayout(layoutBot)

        dockBot.setWidget(self.dataWidget)
        self.addDockWidget(Qt.LeftDockWidgetArea, dockBot)



    def createCentralWidget(self):
        
        self.centralWidget = QWidget()
        
        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        self.toolbar = NavigationToolbar(self.canvas, self)

        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.centralWidget.setLayout(layout)

        self.setCentralWidget(self.centralWidget)

        # test with random data
        self.test()
        
# slots for menu actions

    def export(self):
        filename, filtr = QFileDialog.getSaveFileName(
            self,
            "Choose a file name", '.', "PDF (*.pdf)")
        if not filename:
            return

    def script(self):
        filename, filtr = QFileDialog.getSaveFileName(
            self,
            "Choose a file name", '.', "Python (*.py)")
        if not filename:
            return

    def about(self):
        QMessageBox.about(
            self, "About GetDist GUI",
            "Qt application for GetDist plots.")

# slots for selectWidget 


    def selectRootDir(self):
        
        title = self.tr("Choose an existing case")
        path = os.getcwd()
        dirName = QFileDialog.getExistingDirectory(
            self, title, path,
            QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks)

        if dirName:
            self.rootdir = str(dirName)
            self.lineEditDirectory.setText(self.rootdir)
            
            # 
            fileobj = os.path.join(self.rootdir, "batch.pyobj")
            if os.path.isfile(fileobj):
                self._updateComboBoxParamTag()
            else:
                filesparam = glob.glob(os.path.join(self.rootdir, "*.paramnames"))
                if filesparam:
                    fileparam = os.path.basename(filesparam[0])
                    # read paramnames

    def _updateComboBoxParamTag(self):
        if self.rootdir and os.path.isdir(self.rootdir):
            self.comboBoxParamTag.clear()
            dirs = [ d for d in os.listdir(self.rootdir) if os.path.isdir(os.path.join(self.rootdir, d)) ]
            self.comboBoxParamTag.addItems(dirs)

    def setParamTag(self, strParamTag):
        self.paramTag = str(strParamTag)
        paramDir = os.path.join(self.rootdir, self.paramTag)
        if os.path.isdir(paramDir):
            dirs = [ d for d in os.listdir(paramDir) if os.path.isdir(os.path.join(paramDir, d)) ]
            self.comboBoxDataTag.addItems(dirs)
    
    def setDataTag(self, strDataTag):
        print "strDataTag ", strDataTag
        
    def clickParamX(self, item):
        pass

    def clickParamY(self, item):
        pass

    def statusSelectAllX(self):
        if self.selectAllX.isChecked():
            state = Qt.Checked
        else: 
            state = Qt.Unchecked

        for i in range(self.listParametersX.count()):
            self.listParametersX.item(i).setCheckState(state)

    def statusSelectAllY(self):
        if self.selectAllY.isChecked():
            state = Qt.Checked
        else: 
            state = Qt.Unchecked

        for i in range(self.listParametersY.count()):
            self.listParametersY.item(i).setCheckState(state)


    def plotData(self):
        g = GetDistPlots.GetDistPlotter('PLA/plot_data')
        



   
# slots for selectWidget 

    def openNewFile(self):
        
        title = self.tr("Choose an existing file")
        path = os.getcwd()
        fileName = QFileDialog.getOpenFileName(self, title, path)

        if fileName:
            textFileHandle = open(fileName)
            textFileLines = textFileHandle.read()
            textFileHandle.close()
            
            textBrowser = QTextBrowser(self.tabWidget)
            textBrowser.setPlainText(textFileLines)
            textBrowser.setReadOnly(True)
            
            self.tabWidget.addTab(textBrowser, os.path.basename(str(fileName)))
            self.tabWidget.setTabToolTip(self.tabWidget.currentIndex(), os.path.basename(str(fileName)))

    def closeTab(self, index):
        self.tabWidget.removeTab(index)
                         

# test functions

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
