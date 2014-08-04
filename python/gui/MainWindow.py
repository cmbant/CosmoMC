# -*- coding: utf-8 -*-

import os
import sys
import glob
import signal
import random
import logging

try:
    from PyQt4.QtCore import *
    from PyQt4.QtGui  import *
except ImportError:
    print "Can't import PyQt4.QtCore or PyQt4.QtGui modules." 
    sys.exit()

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import matplotlib.pyplot as plt

import Resources

import GetDistPlots
import MCSamples

# ==============================================================================

class MainWindow(QMainWindow):

    def __init__(self):
        """
        Initialize of GUI components.
        """

        super(MainWindow, self).__init__()

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

        # Internal data 
        self.rootdir = ""
        self.mcsamples = None


    def createActions(self):
        self.exportAct = QAction(QIcon(":/images/file_export.png"),
            "&Export as PDF", self,
            statusTip="Export image as PDF", triggered=self.export)

        self.scriptAct = QAction(QIcon(":/images/file_save.png"), 
            "Save script", self,
            statusTip="Export commands to script", triggered=self.script)

        self.exitAct = QAction(QIcon(":/images/application_exit.png"), 
            "E&xit", self, shortcut="Ctrl+Q",
            statusTip="Exit application",
            triggered=qApp.closeAllWindows)
 
        self.aboutAct = QAction(QIcon(":/images/help_about.png"), 
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

        self.pushButtonSelect = QPushButton(QIcon(":/images/folder_open.png"), 
            "", self.selectWidget)
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
        self.connect(self.listParametersX, SIGNAL("itemChanged(QListWidgetItem *)"), self.itemParamXChanged)

        self.listParametersY = QListWidget(self.selectWidget)
        self.listParametersY.clear()
        self.connect(self.listParametersY, SIGNAL("itemChanged(QListWidgetItem *)"), self.itemParamYChanged)

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

        self.newButton = QPushButton(QIcon(":/images/file_open.png"), "", self.dataWidget)
        self.newButton.setToolTip("Select file to display")
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

        # Canvas Widget that displays the `figure`
        self.canvas = FigureCanvas(self.figure)

        # Navigation widget
        self.toolbar = NavigationToolbar(self.canvas, self)

        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.centralWidget.setLayout(layout)

        self.setCentralWidget(self.centralWidget)

        self.test() # test with random data
        

    # slots for menu actions

    def export(self):
        filename, filter = QFileDialog.getSaveFileName(
            self,
            "Choose a file name", '.', "PDF (*.pdf)")
        if not filename:
            return

    def script(self):
        filename, filter = QFileDialog.getSaveFileName(
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
        """
        Slot function called when pushButtonSelect is pressed.
        """
        
        title = self.tr("Choose an existing case")
        path = os.getcwd()
        dirName = QFileDialog.getExistingDirectory(
            self, title, path,
            QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks)

        #fixme: dialog window won't close immediately (on ubuntu 12.04) 

        if dirName:

            # Root directory
            self.rootdir = str(dirName)
            self.lineEditDirectory.setText(self.rootdir)

            # Root file name
            self.root = MCSamples.GetRootFileName(self.rootdir)            
            logging.debug("self.root = %s"%str(self.root))

            # Get list of chain files in root directory
            chainFiles = MCSamples.GetChainFiles(self.root)

            if len(chainFiles)==0:
                logging.debug("No chain files in %s"%self.rootdir)
                self._updateComboBoxParamTag()
            else:
                # Create instance of MCSamples and read chain files
                self.mcsamples = MCSamples.MCSamples(self.root)
                logging.debug("Read chains from %s"%str(self.root))
                logging.debug("Please wait ...")
                self.mcsamples.loadChains(self.root, chainFiles)
                logging.debug("Done reading")

                # File .paramnames
                filesparam = glob.glob(os.path.join(self.rootdir, "*.paramnames"))
                if filesparam: 
                    # Directory contains file .paramnames
                    if hasattr(self.mcsamples, "paramNames"):
                        logging.debug("Reading .paramnames information")
                        paramNames = self.mcsamples.paramNames.list()

                        # Hide combo boxes and fill list
                        self.comboBoxParamTag.hide()
                        self.comboBoxDataTag.hide()
                        self._updateListParametersX(paramNames)
            
                # File batch.pyobj
                fileobj = os.path.join(self.rootdir, "batch.pyobj")
                if os.path.isfile(fileobj):
                    pass
            
        else:
            logging.debug("No directory specified")


    def _updateComboBoxParamTag(self, listOfParams=[]):
        if self.rootdir and os.path.isdir(self.rootdir):
            self.comboBoxParamTag.clear()
            if not listOfParams:
                listOfParams = [ d for d in os.listdir(self.rootdir) if os.path.isdir(os.path.join(self.rootdir, d)) ]
            self.comboBoxParamTag.addItems(listOfParams)

    def setParamTag(self, strParamTag):
        """
        Slot function called on change of comboBoxParamTag.
        """
        self.paramTag = str(strParamTag)
        paramDir = os.path.join(self.rootdir, self.paramTag)
        if os.path.isdir(paramDir):
            dirs = [ d for d in os.listdir(paramDir) if os.path.isdir(os.path.join(paramDir, d)) ]
            self.comboBoxDataTag.addItems(dirs)
    
    def setDataTag(self, strDataTag):
        """
        Slot function called on change of comboBoxDataTag.
        """
        print "strDataTag ", strDataTag
        

    def _updateListParametersX(self, items):
        """
        """
        logging.debug("Fill ListParametersX with %i items"%(len(items)))
        self.listParametersX.clear()
        for item in items:
            listItem = QListWidgetItem()
            listItem.setText(item)
            listItem.setFlags(listItem.flags() | Qt.ItemIsUserCheckable) 
            listItem.setCheckState(Qt.Unchecked)
            self.listParametersX.addItem(listItem)

    def itemParamXChanged(self, item):
        """
        Slot function called when item in listParametersX changes.
        """
        if item.checkState()==Qt.Checked:
            logging.debug("Change of Item %s: now in check state"%str(item.text()))


    def _updateListParametersY(self, items):
        """
        """
        logging.debug("Fill ListParametersY with %i items"%(len(items)))
        self.listParametersY.clear()
        for item in items:
            listItem = QListWidgetItem()
            listItem.setText(item)
            listItem.setFlags(listItem.flags() | Qt.ItemIsUserCheckable) 
            listItem.setCheckState(Qt.Unchecked)
            self.listParametersY.addItem(listItem)

    def itemParamYChanged(self, item):
        """
        Slot function called when item in listParametersY changes.
        """
        if item.checkState()==Qt.Checked:
            logging.debug("Change of Item %s: now in check state"%str(item.text()))


    def statusSelectAllX(self):
        """
        Slot function called when selectAllX is modified.
        """
        if self.selectAllX.isChecked():
            state = Qt.Checked
        else: 
            state = Qt.Unchecked

        for i in range(self.listParametersX.count()):
            self.listParametersX.item(i).setCheckState(state)

    def statusSelectAllY(self):
        """
        Slot function called when selectAllY is modified.
        """
        if self.selectAllY.isChecked():
            state = Qt.Checked
        else: 
            state = Qt.Unchecked

        for i in range(self.listParametersY.count()):
            self.listParametersY.item(i).setCheckState(state)


    def plotData(self):
        """
        Slot function called when pushButtonPlot is pressed.
        """

        # X items 
        items_x = []
        count = self.listParametersX.count()
        for i in range(count):
            item = self.listParametersX.item(i)
            if item.checkState()==Qt.Checked:
                items_x.append(str(item.text()))

        # X items 
        items_y = []
        count = self.listParametersY.count()
        for i in range(count):
            item = self.listParametersY.item(i)
            if item.checkState()==Qt.Checked:
                items_y.append(str(item.text()))


        from PyQt4.QtCore import pyqtRemoveInputHook
        from pdb import set_trace
        pyqtRemoveInputHook()
        set_trace()

        g = GetDistPlots.GetDistPlotter('PLA/plot_data')
        g.settings.setWithSubplotSize(3.000000)

        # 
        self.plotter = GetDistPlots.GetDistPlotter('main/plot_data')


    #
    # slots for selectWidget (bottom widget on the left)
    #

    def openNewFile(self):
        """
        Slot function called when newButton is pressed.
        """
        
        title = self.tr("Choose an existing file")
        path = os.getcwd()
        fileName = QFileDialog.getOpenFileName(self, title, path)

        if fileName:
            # Read content of file
            textFileHandle = open(fileName)
            textFileLines = textFileHandle.read()
            textFileHandle.close()
            
            # Display contentin read-only mode
            textBrowser = QTextBrowser(self.tabWidget)
            textBrowser.setPlainText(textFileLines)
            textBrowser.setReadOnly(True)

            # Connection of signal
            self.tabWidget.addTab(textBrowser, os.path.basename(str(fileName)))
            self.tabWidget.setTabToolTip(self.tabWidget.currentIndex(), os.path.basename(str(fileName)))

        else:
            # Do nothing
            pass


    def closeTab(self, index):
        """
        Slot function called on a tab widget close event.
        """
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


# ==============================================================================

if __name__ == "__main__":

    # Configure the logging
    level = logging.DEBUG
    FORMAT = '%(asctime).19s [%(levelname)s]\t[%(filename)s:%(lineno)d]\t\t%(message)s'
    logging.basicConfig(level=level, format=FORMAT)      

    app = QApplication(sys.argv)
    mainWin = MainWindow()
    mainWin.show()
    sys.exit(app.exec_())

# ==============================================================================


