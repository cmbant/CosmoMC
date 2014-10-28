# -*- coding: utf-8 -*-

import os
import sys
import glob
import signal
import logging

import matplotlib
#matplotlib.use('Qt4Agg')

try:
    from PySide.QtCore import *
    from PySide.QtGui  import *
    os.environ['QT_API'] = 'pyside'
    matplotlib.rcParams['backend.qt4']='PySide'
    try:
        import Resources_pyside
    except ImportError:
        print "Missing Resources_pyside.py: Run script update_resources.sh"
except ImportError:
    print "Can't import PySide modules."
    sys.exit()

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
#import matplotlib.pyplot as plt

import GetDistPlots
import MCSamples

import batchJob

# ==============================================================================

class MainWindow(QMainWindow):

    def __init__(self):
        """
        Initialize of GUI components.
        """
        super(MainWindow, self).__init__()

        # GUI setup
        self.createDockWidgets()
        self.createActions()
        self.createMenus()
        self.createStatusBar()

        self.setAttribute(Qt.WA_DeleteOnClose)
        self.setWindowTitle("GetDist GUI")
        self.resize(800, 600)

        # Allow to shutdown the GUI with Ctrl+C
        signal.signal(signal.SIGINT, signal.SIG_DFL)

        # Path for .ini file
        self.iniFile = ""

        # Path of root directory
        self.rootdirname = None

        # Root name for chain
        self.rootname = None

        # Dict for other roots
        self.other_rootnames = {}

        # Grid chains parameters
        self.is_grid = False
        self.grid_items = {}
        self.paramTag = ""
        self.dataTag = ""
        self.data2chains = {}

        # Plot instance
        self.plotter = None

        # Script
        self.script = ""


    def createActions(self):
        """
        Create Qt actions used in GUI.
        """
        self.exportAct = QAction(QIcon(":/images/file_export.png"),
                                 "&Export as PDF", self,
                                 statusTip="Export image as PDF",
                                 triggered=self.export)

        self.scriptAct = QAction(QIcon(":/images/file_save.png"),
                                 "Save script", self,
                                 statusTip="Export commands to script",
                                 triggered=self.saveScript)

        self.exitAct = QAction(QIcon(":/images/application_exit.png"),
                               "E&xit", self,
                               shortcut="Ctrl+Q",
                               statusTip="Exit application",
                               triggered=self.close)

        self.statsAct = QAction(QIcon(":/images/view_text.png"),
                               "Marge Stats", self,
                               shortcut="",
                               statusTip="Show Marge Stats",
                               triggered=self.showMargeStats)

        self.aboutAct = QAction(QIcon(":/images/help_about.png"),
                                "&About", self,
                                statusTip="Show About box",
                                triggered=self.about)

    def createMenus(self):
        """
        Create Qt menus.
        """
        self.fileMenu = self.menuBar().addMenu("&File")
        self.fileMenu.addAction(self.exportAct)
        self.fileMenu.addAction(self.scriptAct)
        self.separatorAct = self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.exitAct)

        self.menuBar().addSeparator()
        self.dataMenu = self.menuBar().addMenu("&Data")
        self.dataMenu.addAction(self.statsAct)

        self.menuBar().addSeparator()

        self.windowMenu = self.menuBar().addMenu("&Windows")
        self.windowMenu.addAction(self.dockTop.toggleViewAction())

        self.menuBar().addSeparator()

        self.helpMenu = self.menuBar().addMenu("&Help")
        self.helpMenu.addAction(self.aboutAct)

    def createStatusBar(self):
        """
        Create Qt status bar.
        """
        self.statusBar().showMessage("Ready", 2000)

    def createDockWidgets(self):
        """
        Create dockable widget and its content.
        """
        self.dockTop = QDockWidget(self.tr(""), self)
        self.dockTop.setAllowedAreas(Qt.LeftDockWidgetArea)

        self.selectWidget = QWidget(self.dockTop)

        self.lineEditDirectory = QLineEdit(self.selectWidget)
        self.lineEditDirectory.clear()
        self.lineEditDirectory.setReadOnly(True)

        self.pushButtonSelect = QPushButton(QIcon(":/images/file_add.png"),
                                            "", self.selectWidget)
        self.pushButtonSelect.setToolTip("Choose root directory")
        self.connect(self.pushButtonSelect, SIGNAL("clicked()"),
                     self.selectRootDirName)
        shortcut = QShortcut(QKeySequence(self.tr("Ctrl+O")), self)
        self.connect(shortcut, SIGNAL("activated()"), self.selectRootDirName)

        self.listRoots = QListWidget(self.selectWidget)
        self.listRoots.setMaximumSize(QSize(16777215, 120))
        self.listRoots.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.listRoots.setSelectionMode(QAbstractItemView.SingleSelection)
        self.connect(self.listRoots,
                     SIGNAL("itemChanged(QListWidgetItem *)"),
                     self.updateListRoots)

        self.pushButtonRemove = QPushButton(QIcon(":/images/file_remove.png"),
                                            "", self.selectWidget)
        self.pushButtonRemove.setToolTip("Remove a root directory")
        self.connect(self.pushButtonRemove, SIGNAL("clicked()"),
                     self.removeOtherRoot)


        self.comboBoxParamTag = QComboBox(self.selectWidget)
        self.comboBoxParamTag.clear()
        self.connect(self.comboBoxParamTag,
                     SIGNAL("activated(const QString&)"), self.setParamTag)

        self.comboBoxDataTag = QComboBox(self.selectWidget)
        self.comboBoxDataTag.clear()
        self.connect(self.comboBoxDataTag,
                     SIGNAL("activated(const QString&)"), self.setDataTag)

        self.listParametersX = QListWidget(self.selectWidget)
        self.listParametersX.clear()
        self.connect(self.listParametersX,
                     SIGNAL("itemChanged(QListWidgetItem *)"),
                     self.itemParamXChanged)

        self.listParametersY = QListWidget(self.selectWidget)
        self.listParametersY.clear()
        self.connect(self.listParametersY,
                     SIGNAL("itemChanged(QListWidgetItem *)"),
                     self.itemParamYChanged)

        self.selectAllX = QCheckBox("Select All", self.selectWidget)
        self.selectAllX.setCheckState(Qt.Unchecked)
        self.connect(self.selectAllX, SIGNAL("clicked()"),
                     self.statusSelectAllX)

        self.selectAllY = QCheckBox("Select All", self.selectWidget)
        self.selectAllY.setCheckState(Qt.Unchecked)
        self.connect(self.selectAllY, SIGNAL("clicked()"),
                     self.statusSelectAllY)

        self.toggleFilled = QRadioButton("Filled")
        self.toggleLine = QRadioButton("Line")
        self.toggleColor = QRadioButton("Color by:")
        self.comboBoxColor = QComboBox(self)
        self.comboBoxColor.clear()
        self.toggleFilled.setChecked(True)

        self.trianglePlot = QCheckBox("Triangle Plot", self.selectWidget)
        self.trianglePlot.setCheckState(Qt.Unchecked)

        self.pushButtonPlot = QPushButton("Make plot", self.selectWidget)
        self.connect(self.pushButtonPlot, SIGNAL("clicked()"), self.plotData)

        # Graphic Layout
        layoutTop = QGridLayout()
        layoutTop.setSpacing(5)
        layoutTop.addWidget(self.lineEditDirectory, 0, 0, 1, 3)
        layoutTop.addWidget(self.pushButtonSelect,  0, 3, 1, 1)
        layoutTop.addWidget(self.listRoots,         1, 0, 2, 3)
        layoutTop.addWidget(self.pushButtonRemove,  1, 3, 1, 1)
        layoutTop.addWidget(self.comboBoxParamTag,  3, 0, 1, 2)
        layoutTop.addWidget(self.comboBoxDataTag,   3, 2, 1, 2)
        layoutTop.addWidget(self.selectAllX,        4, 0, 1, 2)
        layoutTop.addWidget(self.selectAllY,        4, 2, 1, 2)
        layoutTop.addWidget(self.listParametersX,   5, 0, 5, 2)
        layoutTop.addWidget(self.listParametersY,   5, 2, 1, 2)
        layoutTop.addWidget(self.toggleFilled,      6, 2, 1, 2)
        layoutTop.addWidget(self.toggleLine,        7, 2, 1, 2)
        layoutTop.addWidget(self.toggleColor,       8, 2, 1, 1)
        layoutTop.addWidget(self.comboBoxColor,     8, 3, 1, 1)
        layoutTop.addWidget(self.trianglePlot,      9, 2, 1, 2)
        layoutTop.addWidget(self.pushButtonPlot,   11, 0, 1, 4)
        self.selectWidget.setLayout(layoutTop)

        self.listRoots.hide()
        self.pushButtonRemove.hide()
        self.comboBoxParamTag.hide()
        self.comboBoxDataTag.hide()

        # Minimum size for initial resize of dock widget
        #self.selectWidget.setMinimumSize(250, 250)
        self.selectWidget.setMaximumSize(600, 1200)

        self.dockTop.setWidget(self.selectWidget)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.dockTop)

    def setIniFile(self, iniFile):
        logging.debug("ini file is %s"%iniFile)
        self.iniFile = iniFile

    # slots for menu actions

    def export(self):
        """
        Callback for action 'Export as PDF'.
        """
        if self.plotter:
            filename, filt = QFileDialog.getSaveFileName(
                self, "Choose a file name", '.', "PDF (*.pdf)")
            if not filename: return
            filename = str(filename)

            logging.debug("Saving PDF in file %s"%filename)
            self.plotter.export(filename)
        else:
            #logging.warning("No plotter data to export")
            QMessageBox.warning(self, "Export", "No plotter data to export")

    def saveScript(self):
        """
        Callback for action 'Save script'.
        """
        if self.script=='':
            QMessageBox.warning(self, "Script", "No script to save")
            #logging.warning("Script is empty!")
            return

        filename, filt = QFileDialog.getSaveFileName(
            self, "Choose a file name", '.', "rPython (*.py)")
        if not filename: return
        filename = str(filename)
        logging.debug("Export script to %s"%filename)
        f = open(filename, 'w')
        f.write(self.script)
        f.close()

    def showMargeStats(self):
        """
        Callback for action 'Show Marge Stats'.
        """
        if self.rootname is not None:
            rootname = self.rootname
        else:
            QMessageBox.warning(self, "Marge Stats", "No rootname. Can't show marge stats")
            #logging.warning("No rootname. Can't show marge stats")
            return

        if self.plotter is None:
            QMessageBox.warning(self, "Marge Stats", "No plotter. Can't show marge stats")
            #logging.warning("No plotter. Can't show marge stats")
            self.statusBar().showMessage("No data available", 2000)
            return

        text = self.plotter.sampleAnalyser.getMargeStats(rootname)
        text = text.replace('\t', '\t\t') # uses double tab
        dlg = DialogMargeStats(self, text)
        dlg.exec_()

    def about(self):
        """
        Callback for action 'About'.
        """
        QMessageBox.about(
            self, "About GetDist GUI",
            "Qt application for GetDist plots.")


    # slots for selectWidget

    def selectRootDirName(self):
        """
        Slot function called when pushButtonSelect is pressed.
        """

        # If a rootdir is defined, add another root
        if self.rootdirname is not None and self.is_grid==False:
            self.addRoot()
            return

	settings = QSettings('cosmomc', 'gui')
        last_dir = settings.value('lastSearchDirectory')
        last_dir = ''

	# Search in current directory, if no previous path available
        if not last_dir: last_dir = os.getcwd()

        title = self.tr("Choose an existing case")
        dirName = QFileDialog.getExistingDirectory(
            self, title, last_dir,
            QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks)
        dirName = str(dirName)
        logging.debug("dirName: %s"%dirName)

        # No directory selected
        if dirName is None or dirName=='':
            return

        settings.setValue('lastSearchDirectory', dirName)

        # Root directory
        self.rootdirname = dirName
        self.lineEditDirectory.setText(self.rootdirname)

        # Grid chains
        filebatch = os.path.join(dirName, "batch.pyobj")
        if os.path.isfile(filebatch):
            self._readGridChain(self.rootdirname)
            return

        # Root file name
        rootname = self._getRootFileName()
        if rootname=='':
            QMessageBox.warning(self, "Select root", "Root file name is not defined")
            #logging.warning("Root file name is not defined")
            return
        else:
            self.rootname = rootname
            logging.info("rootname: %s"%self.rootname)

        # Get chain files
        chainFiles = MCSamples.GetChainFiles(self.rootname)
        if len(chainFiles)==0:
            logging.debug("No chain files in %s"%self.rootdirname)
            self._updateComboBoxParamTag()
            self._updateComboBoxColor()
            return

        self.plotter = GetDistPlots.GetDistPlotter(mcsamples=True,
                                                   ini_file=self.iniFile)
        logging.debug("Read chains in %s"%str(self.rootname))
        self.plotter.sampleAnalyser.addRoot(self.rootname)
        self._updateParameters()

        # Hide combo boxes and fill list
        self.comboBoxParamTag.hide()
        self.comboBoxDataTag.hide()
        self.listRoots.show()
        self.pushButtonRemove.show()



    def _updateParameters(self):
        paramSets = []

        if self.rootname is not None:
            params = self.plotter.sampleAnalyser.usedParamsForRoot(self.rootname)
            paramSets.append(set(params))

        for k, values in self.other_rootnames.items():
            rootname, state = values
            if state:
                params = self.plotter.sampleAnalyser.usedParamsForRoot(rootname)
                paramSets.append(set(params))

        if paramSets:
            intersect = set.intersection(*paramSets)
            paramNames = list(intersect)
            self._updateListParametersX(paramNames)
            self._updateListParametersY(paramNames)
            self._updateComboBoxColor(paramNames)
        else:
            logging.warning("No parameters found")


    def _readGridChain(self, batchPath):
        """
        Setup of a grid chain.
        """
        self.is_grid = True
        logging.debug("Read grid chain in %s"%batchPath)
        batch = batchJob.readobject(batchPath)
        items = dict()
        for jobItem in batch.items(True, True):
            if jobItem.chainExists():
                if not jobItem.paramtag in items: items[jobItem.paramtag] = []
                items[jobItem.paramtag].append((jobItem.datatag, jobItem.chainRoot))
        logging.debug("Found %i names for grid"%len(items.keys()))
        self.grid_items = items
        self.comboBoxParamTag.show()
        self.comboBoxDataTag.show()
        self.listRoots.show()
        self.pushButtonRemove.show()
        self._updateComboBoxParamTag(self.grid_items.keys())

    def _getFileParam(self, paramTag=None, dataTag=None):
        """
        Get corresponding .paramnames file.
        """
        fileparam = ""
        if self.rootdirname is not None:
            if paramTag is not None and dataTag is not None:
                rootdirname = os.path.join(self.rootdirname, paramTag, dataTag)
            else:
                rootdirname = self.rootdirname
            filesparam = MCSamples.GetParamNamesFiles(rootdirname)
            if len(filesparam)>1:
                d = dict()
                for fileparam in filesparam:
                    d[os.path.basename(fileparam)] = fileparam
                fn, ok = QInputDialog.getItem(
                    self, "Select file", "Param name:", d.keys(), 0, False)
                if ok and fn:
                    fileparam = d[fn]
            elif len(filesparam)==1:
                fileparam = filesparam[0]
        return fileparam

    def _getRootFileName(self, paramTag=None, dataTag=None):
        """
        Get the root file name.
        """
        rootname = ''
        fileparam = self._getFileParam(paramTag, dataTag)
        if fileparam<>'':
            rootname = str(fileparam).replace(".paramnames", "")
        else:
            if paramTag is not None and dataTag is not None:
                rootdirname = os.path.join(self.rootdirname, paramTag, dataTag)
            else:
                rootdirname = self.rootdirname
            rootname = MCSamples.GetRootFileName(rootdirname)
        return rootname

    def addRoot(self):
        """
        Add another root file name.
        """
	settings = QSettings('cosmomc', 'gui')
        last_dir = settings.value('lastSearchDirectory')

	# Search in current directory, if no previous path available
        if not last_dir: last_dir = os.getcwd()

        filt = "Param names (*.paramnames)"
        title = self.tr("Choose a chain root")
        fileName, filt = QFileDialog.getOpenFileName(
            self, title, last_dir, filt)
        fileName = str(fileName)

        if fileName:

            dirName = os.path.dirname(fileName)
            settings.setValue('lastSearchDirectory', dirName)

            root = fileName.replace('.paramnames', '')
            basename = os.path.basename(root)
            if root==self.rootname:
                logging.warning("Same root as main root")
                return
            elif not self.other_rootnames.has_key(basename):
                if self.plotter is not None:
                    self.plotter.sampleAnalyser.addRoot(root)
                    logging.debug("Add root %s"%basename)
                    self.other_rootnames[basename] = [root, True]
                    item = QListWidgetItem(self.listRoots)
                    item.setText(basename)
                    item.setCheckState(Qt.Checked)
                    idx = self.listRoots.count()
                    self.listRoots.insertItem(idx, item)
                    self._updateParameters()
                else:
                    logging.warning("No plotter instance")

    def updateListRoots(self, item):
        logging.debug("updateListRoots")
        nitems = self.listRoots.count()
        for i in range(nitems):
            item = self.listRoots.item(i)
            root = item.text()
            state = item.checkState()==Qt.Checked
            if self.other_rootnames.has_key(root):
                self.other_rootnames[root][1] = state
            else:
                logging.warning('No root found for %s'%root)
        self._updateParameters()

    def removeOtherRoot(self):
        logging.debug("Remove root")
        for i in range(self.listRoots.count()):
            item = self.listRoots.item(i)
            if item.isSelected():
                root = str(item.text())
                logging.debug("Remove root %s"%root)
                self.plotter.sampleAnalyser.removeOtherRoot(root)
                self.listRoots.takeItem(i)
                if self.other_rootnames.has_key(root):
                    del self.other_rootnames[root]

    def getOtherRoots(self):
        logging.debug("Get status for other roots")
        for i in range(self.listRoots.count()):
            item = self.listRoots.item(i)
            root = str(item.text())
            state = (item.checkState()==Qt.Checked)
            if self.other_rootnames.has_key(root):
                self.other_rootnames[root][1] = state
                logging.debug("status for '%s': %s"%(root, state))
            else:
                logging.warning("No key for root '%s'"%root)

    def _updateComboBoxParamTag(self, listOfParams=[]):
        if self.rootdirname and os.path.isdir(self.rootdirname):
            self.comboBoxParamTag.show()
            self.comboBoxParamTag.clear()
            if not listOfParams:
                listOfParams = [
                    d for d in os.listdir(self.rootdirname)
                    if os.path.isdir(os.path.join(self.rootdirname, d)) ]
            self.comboBoxParamTag.addItems(listOfParams)
            self.setParamTag(str(self.comboBoxParamTag.itemText(0)))

    def setParamTag(self, strParamTag):
        """
        Slot function called on change of comboBoxParamTag.
        """
        self.paramTag = str(strParamTag)
        logging.debug("Param: %s"%self.paramTag)
        if self.is_grid:
            if self.grid_items.has_key(self.paramTag):
                values = self.grid_items[self.paramTag]
                datas = [ (v[0], v[1]) for v in values ]
                self._updateComboBoxDataTag(datas)
            else:
                logging.warning("Grid defined but no data found for param %s"%self.paramTag)
        else:
            logging.warning("No grid defined and no data found for param %s"%self.paramTag)
            dir_param = os.path.join(self.rootdirname, self.paramTag)
            if os.path.isdir(dir_param):
                subdirs = [ d for d in os.listdir(dir_param)
                            if os.path.isdir(os.path.join(dir_param, d)) ]
                self._updateComboBoxDataTag(subdirs)
            else:
                logging.warning("No param found")

    def _updateComboBoxDataTag(self, listOfDatas=[]):
        logging.debug("Update data with %i values"%len(listOfDatas))
        if self.rootdirname and os.path.isdir(self.rootdirname):
            self.comboBoxDataTag.show()
            self.comboBoxDataTag.clear()
            if listOfDatas:
                self.data2chains = {}
                for data in listOfDatas:
                    if type(data)==type(()): # (dataTag, chainRoot)
                        #self.comboBoxDataTag.addItem(data[0], data[1])
                        self.comboBoxDataTag.addItem(data[0])
                        self.data2chains[data[0]] = data[1]
                        logging.debug("Insert %s (%s)"%(str(data[0]), str(data[1])))
                    else:
                        self.comboBoxDataTag.addItem(data) # dataTag
            else:
                logging.warning("No datas to display")

    def setDataTag(self, strDataTag):
        """
        Slot function called on change of comboBoxDataTag.
        """
        self.dataTag = str(strDataTag)
        logging.debug("Data: %s"%strDataTag)

        if self.data2chains.has_key(strDataTag):
            rootname = self.data2chains.get(strDataTag)
            logging.debug("rootname: %s"%rootname)
        else:
            logging.warning("Root file name is not defined")
            return

        basename = os.path.basename(rootname)
        if self.other_rootnames.has_key(basename):
            logging.warning("Root name '%s' allready selected"%basename)
            return

        if self.plotter is None:
            self.plotter = GetDistPlots.GetDistPlotter(mcsamples=True,
                                                       ini_file=self.iniFile)
        self.statusBar().showMessage("Read chains in %s"%str(rootname))
        self.plotter.sampleAnalyser.addRoot(rootname)
        self.other_rootnames[basename] = [rootname, True]
        item = QListWidgetItem(self.listRoots)
        item.setText(basename)
        item.setCheckState(Qt.Checked)
        idx = self.listRoots.count()
        self.listRoots.insertItem(idx, item)
        self._updateParameters()


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

    def _updateComboBoxColor(self, listOfParams=[]):
        if self.rootdirname and os.path.isdir(self.rootdirname):
            self.comboBoxColor.clear()
            if not listOfParams:
                listOfParams = [
                    d for d in os.listdir(self.rootdirname)
                    if os.path.isdir(os.path.join(self.rootdirname, d)) ]
            self.comboBoxColor.addItems(listOfParams)


    def plotData(self):
        """
        Slot function called when pushButtonPlot is pressed.
        """
        # Ensure at least 1 root name specified
        self.getOtherRoots()
        if self.rootname is None:
            if len(self.other_rootnames)==0:
                QMessageBox.warning(self, "Plot data", "No root defined")
                #logging.warning("No rootname defined")
                return
            else:
                nroots = 0
                for basename, values in self.other_rootnames.items():
                    other_rootname, state = values
                    if state: nroots += 1
                if nroots==0:
                     logging.warning("No rootname selected")
                     QMessageBox.warning(self, "Plot data", "No root selected")
                     return

        if self.plotter is None:
            QMessageBox.warning(self, "Plot data", "No GetDistPlotter instance")
            #logging.warning("No GetDistPlotter instance")
            return

        self.plotter.settings.setWithSubplotSize(3.000000)

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


        # Script
        self.script  = ""
        self.script += "import GetDistPlots, os\n"
        self.script += "g=GetDistPlots.GetDistPlotter(mcsamples=True, ini_file='%s')\n"%self.iniFile
        self.script += "g.settings.setWithSubplotSize(3.000000)\n"
        self.script += "outdir='.'\n"
        self.script += "dict_roots={}\n"

        roots = []
        # Main root name
        if self.rootname is not None:
            basename = os.path.basename(self.rootname)
            roots.append(basename)
            self.script += "dict_roots['%s'] = '%s'\n"%(basename, self.rootname)

        # Other root names
        for basename, values in self.other_rootnames.items():
            other_rootname, state = values
            if state:
                roots.append(os.path.basename(other_rootname))
                self.script += "dict_roots['%s'] = '%s'\n"%(basename, other_rootname)
            else:
                logging.debug("root '%s' not selected"%basename)

        self.script += "g.sampleAnalyser.addRoots(dict_roots.values())\n"
        self.script += "roots=dict_roots.keys()\n"

        logging.debug("Plotting with roots = %s"%str(roots))

        if self.trianglePlot.isChecked():
            if len(items_x)>1:
                params = items_x
                logging.debug("Triangle plot ")
                logging.debug("params = %s"%str(params))
                self.script += "params = %s\n"%str(params)
                try:
                    self.plotter.triangle_plot(roots, params)
                except:
                    QMessageBox.critical(
                        self, "Triangle plot",
                        "Error for command 'triangle_plot(roots, params)' with roots=%s \nparams=%s \n"%(str(roots), str(params)))
                    return

                self.updatePlot()
                self.script += "g.triangle_plot(roots, params)\n"
            else:
                QMessageBox.warning(self, "Triangle plot", "Need more than 1 x-param")
                #self.statusBar().showMessage("Need more than 1 x-param", 2000)

        elif len(items_x)>0 and len(items_y)==0:
            params = items_x
            logging.debug("1D plot")
            logging.debug("params = %s"%str(params))
            try:
                self.plotter.plots_1d(roots, params=params)
            except:
                QMessageBox.critical(
                        self, "Plot 1D",
                        "Error for command 'plots_1d(roots, params=params)' with roots=%s \nparams=%s \n"%(str(roots), str(params)))
                return
            self.updatePlot()
            self.script += "params=%s\n"%str(params)
            self.script += "g.plots_1d(roots, params=params)\n"

        elif len(items_x)>0 and len(items_y)>0:
            # only one X item selected
            if self.toggleFilled.isChecked() or self.toggleLine.isChecked():
                filled=self.toggleFilled.isChecked()
                logging.debug("2D plot ")
                #
                pairs = []
                self.script += "pairs = []\n"
                item_x = items_x[0]
                for item_y in items_y:
                    pairs.append([item_x, item_y])
                    self.script += "pairs.append(['%s','%s'])\n"%(item_x, item_y)

                logging.debug("pairs  = %s"%str(pairs))
                self.script += "filled=%s\n"%filled
                try:
                    self.plotter.plots_2d(roots, param_pairs=pairs, filled=filled)
                except:
                    QMessageBox.critical(
                        self, "Plot 2D",
                        "Error for command 'plots_2d(roots, param_pairs=pairs, filled=filled)' with roots=%s \npairs=%s \n"%(str(roots), str(pairs)))
                    return
                self.updatePlot()
                self.script += "g.plots_2d(oots, param_pairs=pairs, filled=filled)\n"

            if self.toggleColor.isChecked():
                sets = []
                sets.append(items_x)
                x = items_x[0]
                y = items_y[0]
                color = str(self.comboBoxColor.currentText())
                logging.debug("3D plot")
                sets = [[x, y, color]]
                logging.debug("sets = %s"%str(sets))
                self.script += "sets = []\n"
                self.script += "sets.append(['%s', '%s', '%s'])\n"%(x, y, color)
                try:
                    self.plotter.plots_3d(roots, sets)
                except:
                    QMessageBox.critical(
                        self, "Plot 3D",
                        "Error for command 'plots_3d(roots, sets)' with roots=%s \nsets=%s \n"%(str(roots), str(sets)))
                    return

                self.updatePlot()
                self.script += "g.plots_3d(roots, sets)\n"

        else:
            text  = ""
            text += "Wrong parameters. Usage:\n"
            text += "\n"
            text += "Tri: box checked + X parameters\n"
            text += "1D : X parameter(s)\n"
            text += "2D : X parameter + Y parameter(s) + Filled or Line\n"
            text += "3D : X parameter + Y parameter + Color parameter\n"
            text += "\n"

            QMessageBox.warning(self, "Plot", text)
            self.script = ""
            return

        if self.rootname is not None:
            basename = os.path.basename(self.rootname)
            self.script += "g.export(os.path.join(outdir,'%s.pdf'))\n"%basename

    def updatePlot(self):

        if self.plotter.fig is None:
            logging.debug("Define an empty central widget")
            self.centralWidget = QWidget()
            self.setCentralWidget(self.centralWidget)

            #self.figure = None
            self.canvas = None
        else:
            logging.debug("Define new canvas in central widget")
            self.centralWidget = QWidget()
            self.setCentralWidget(self.centralWidget)

            self.canvas = FigureCanvas(self.plotter.fig)
            self.toolbar = NavigationToolbar(self.canvas, self)
            layout = QVBoxLayout()
            layout.addWidget(self.toolbar)
            layout.addWidget(self.canvas)
            self.centralWidget.setLayout(layout)

            self.canvas.draw()

# ==============================================================================

class DialogMargeStats(QDialog):

    def __init__(self, parent=None, text=""):
        QDialog.__init__(self, parent)

        self.textBrowser = QTextEdit(self)

        layout = QGridLayout()
        layout.addWidget(self.textBrowser, 0, 0)
        self.setLayout(layout)

        self.setWindowTitle(self.tr("Dialog for MargeStats"))

        if (text):
            self.textBrowser.setPlainText(text)
            self.textBrowser.setReadOnly(True)

# ==============================================================================

if __name__ == "__main__":

    app = QApplication(sys.argv)
    mainWin = MainWindow()
    mainWin.show()
    sys.exit(app.exec_())

# ==============================================================================


