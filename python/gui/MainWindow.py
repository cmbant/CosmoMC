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
    os.environ['QT_API'] = 'pyqt'
    try:
        import Resources_pyqt
    except ImportError:
        print "Missing file Resources_pyqt.py: Run script update_resources.sh"
except ImportError:
    try:
        from PySide.QtCore import *
        from PySide.QtGui  import *
        os.environ['QT_API'] = 'pyside'
        try:
            import Resources_pyside
        except ImportError:
            print "Missing file Resources_pyside.py: Run script update_resources.sh"
    except ImportError:
        print "Can't import PyQt4 or PySide modules." 
        sys.exit()

import matplotlib
matplotlib.use('Qt4Agg')

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import matplotlib.pyplot as plt


import GetDistPlots
import MCSamples

# ==============================================================================

class MainWindow(QMainWindow):

    def __init__(self):
        """
        Initialize of GUI components.
        """

        # Configure the logging
        level = logging.DEBUG
        FORMAT = '%(asctime).19s [%(levelname)s]\t[%(filename)s:%(lineno)d]\t\t%(message)s'
        logging.basicConfig(level=level, format=FORMAT)      


        super(MainWindow, self).__init__()

        self.setAttribute(Qt.WA_DeleteOnClose)
 
        self.createDockWidgets()
        self.createActions()
        self.createMenus()
        #self.createToolBars()
        self.createStatusBar()
        self.createCentralWidget()

        self.setWindowTitle("GetDist GUI")
        self.resize(800, 600)

        # Allow to shutdown the GUI with Ctrl+C
        signal.signal(signal.SIGINT, signal.SIG_DFL)

        # Root directory
        self.rootdir = "" 

        self.plotter = None


    def createActions(self):
        self.exportAct = QAction(QIcon(":/images/file_export.png"),
                                 "&Export as PDF", self,
                                 statusTip="Export image as PDF", 
                                 triggered=self.export)

        self.scriptAct = QAction(QIcon(":/images/file_save.png"), 
                                 "Save script", self,
                                 statusTip="Export commands to script", 
                                 triggered=self.script)

        self.exitAct = QAction(QIcon(":/images/application_exit.png"), 
                               "E&xit", self, 
                               shortcut="Ctrl+Q",
                               statusTip="Exit application",
                               triggered=self.close)
 
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
 
        self.windowMenu = self.menuBar().addMenu("&Windows")
        self.windowMenu.addAction(self.dockTop.toggleViewAction())
        self.windowMenu.addAction(self.dockBot.toggleViewAction())
                                  
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
        self.dockTop = QDockWidget(self.tr("Filters"), self)
        self.dockTop.setAllowedAreas(Qt.LeftDockWidgetArea)

        self.selectWidget = QWidget(self.dockTop)

        self.lineEditDirectory = QLineEdit(self.selectWidget)
        self.lineEditDirectory.clear()
        self.lineEditDirectory.setReadOnly(True)

        self.pushButtonSelect = QPushButton(QIcon(":/images/folder_open.png"), 
            "", self.selectWidget)
        self.pushButtonSelect.setToolTip("Choose root directory")
        self.connect(self.pushButtonSelect, SIGNAL("clicked()"), self.selectRootDir)
        shortcut = QShortcut(QKeySequence(self.tr("Ctrl+O")), self)
        self.connect(shortcut, SIGNAL("activated()"), self.selectRootDir)

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

        self.toggleFilled = QRadioButton("Filled")
        self.toggleLine = QRadioButton("Line")
        self.toggleColor = QRadioButton("Color by:")
        self.lineEditColor = QLineEdit() 
        # color values ???
        self.toggleFilled.setChecked(True)
        

        self.trianglePlot = QCheckBox("Triangle Plot", self.selectWidget)
        self.trianglePlot.setCheckState(Qt.Unchecked)
        #self.connect(self.trianglePlot, SIGNAL("clicked()"), self.setTrianglePlot)

        self.pushButtonPlot = QPushButton("Make plot", self.selectWidget)
        self.connect(self.pushButtonPlot, SIGNAL("clicked()"), self.plotData)

        # Graphic Layout
        layoutTop = QGridLayout()
        layoutTop.setSpacing(5)
        layoutTop.addWidget(self.lineEditDirectory, 0, 0, 1, 3)
        layoutTop.addWidget(self.pushButtonSelect,  0, 3, 1, 1)
        layoutTop.addWidget(self.comboBoxParamTag,  1, 0, 1, 2)
        layoutTop.addWidget(self.comboBoxDataTag,   1, 2, 1, 2)
        layoutTop.addWidget(self.selectAllX,        2, 0, 1, 2)
        layoutTop.addWidget(self.selectAllY,        2, 2, 1, 2)
        layoutTop.addWidget(self.listParametersX,   3, 0, 5, 2)
        layoutTop.addWidget(self.listParametersY,   3, 2, 1, 2)
        layoutTop.addWidget(self.toggleFilled,      4, 2, 1, 2)
        layoutTop.addWidget(self.toggleLine,        5, 2, 1, 2)
        layoutTop.addWidget(self.toggleColor,       6, 2, 1, 1)
        layoutTop.addWidget(self.lineEditColor,     6, 3, 1, 1)
        layoutTop.addWidget(self.trianglePlot,      7, 2, 1, 2)
        layoutTop.addWidget(self.pushButtonPlot,    9, 0, 1, 4)
        self.selectWidget.setLayout(layoutTop)
        
        # Minimum size for initial resize of dock widget
        #self.selectWidget.setMinimumSize(250, 250)
        self.selectWidget.setMaximumSize(350, 350)

        self.dockTop.setWidget(self.selectWidget)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.dockTop)

        # Bottom: Tab Widget  
        self.dockBot = QDockWidget(self.tr("Data Files"), self)
        self.dockBot.setAllowedAreas(Qt.LeftDockWidgetArea)

        self.dataWidget = QWidget(self.dockBot)

        self.newButton = QPushButton(QIcon(":/images/file_open.png"), "", self.dataWidget)
        self.newButton.setToolTip("Select file to display")
        self.connect(self.newButton, SIGNAL("clicked()"), self.openNewFile)

        spacer = QSpacerItem(10, 10, QSizePolicy.Expanding, QSizePolicy.Ignored)

        self.tabWidget = QTabWidget(self.dockBot)
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

        self.dockBot.setWidget(self.dataWidget)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.dockBot)

    def createCentralWidget(self):
        
        self.centralWidget = QWidget()
        
        # a figure instance to plot on
        self.figure = plt.figure()

        # Canvas Widget that displays the 'figure'
        self.canvas = FigureCanvas(self.figure)

        # Navigation widget
        self.toolbar = NavigationToolbar(self.canvas, self)

        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.centralWidget.setLayout(layout)

        self.setCentralWidget(self.centralWidget)

        #self.test() # test with random data
        

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
            chainFiles = MCSamples.GetChainFiles(self.root)
            if len(chainFiles)==0:
                logging.debug("No chain files in %s"%self.rootdir)
                self._updateComboBoxParamTag()
            else:

                self.plotter = GetDistPlots.GetDistPlotter(chain_dir=self.rootdir)
                self.statusBar().showMessage("Reading chain files in %s"%str(self.root))

                # File .paramnames
                filesparam = glob.glob(os.path.join(self.rootdir, "*.paramnames"))
                if filesparam: 
                    # Directory contains file .paramnames
                    paramNames = self.plotter.sampleAnalyser.paramsForRoot(self.rootdir).list()

                    # Hide combo boxes and fill list
                    self.comboBoxParamTag.hide()
                    self.comboBoxDataTag.hide()
                    self._updateListParametersX(paramNames)
            
                # File batch.pyobj
                fileobj = os.path.join(self.rootdir, "batch.pyobj")
                if os.path.isfile(fileobj):
                    pass

                self.statusBar().showMessage("Ready", 2000)
           
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


#    def setTrianglePlot(self):
#        if self.trianglePlot.isChecked():
#            pass


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
            
        
        logging.debug("Create GetDistPlotter with arg %s"%self.rootdir)
        self.plotter.settings.setWithSubplotSize(3.000000)
        roots = [os.path.basename(self.root)]

        #
        if self.trianglePlot.isChecked(): 
            if len(items_x)>1:
                params = items_x
                logging.debug("Triangle plot ")
                logging.debug("roots  = %s"%str(roots))
                logging.debug("params = %s"%str(params))
                self.plotter.triangle_plot(roots, params)
            else:
                self.statusBar().showMessage("Need more than 1 x-param", 2000)
                
        elif len(items_x)>0 and len(items_y)==0:
            params = items_x
            logging.debug("1D plot ")
            logging.debug("roots = %s"%str(roots))
            logging.debug("params = %s"%str(params))
            self.plotter.plots_1d(roots, params=params)
            self.updatePlot()
            
        elif len(items_x)>0 and len(items_y)>0:

            if self.toggleFilled.isChecked() or self.toggleLine.isChecked():
                logging.debug("2D plot ")
                logging.debug("roots  = %s"%str(roots))
                logging.debug("param1 = %s"%str(items_x))
                logging.debug("param2 = %s"%str(items_y))
                self.plotter.plots_2d(roots, param1=items_x, param2=items_y)

            if self.toggleColor.isChecked():
                sets = []
                sets.append(items_x)
                x = items_x[0]
                y = items_y[0]
                color = str(self.lineEditColor.displayText())
                logging.debug("3D plot ")
                logging.debug("roots = %s"%str(roots))
                self.plotter.plot_3d(roots, [x, y, color])

        logging.debug("End of function")
              
    
    def updatePlot(self):
        logging.debug("Update plot ") 
        ax = self.figure.add_subplot(111)
        ax.set_figure(self.plotter.fig)
        self.canvas.draw()


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

    app = QApplication(sys.argv)
    mainWin = MainWindow()
    mainWin.show()
    sys.exit(app.exec_())

# ==============================================================================


