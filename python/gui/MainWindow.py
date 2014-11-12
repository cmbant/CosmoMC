# -*- coding: utf-8 -*-

import os
import sys
#import glob
import signal
import logging

import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

try:
    from PySide.QtCore import *
    from PySide.QtGui  import *
    os.environ['QT_API'] = 'pyside'
    matplotlib.rcParams['backend.qt4'] = 'PySide'
    try:
        import Resources_pyside
    except ImportError:
        print "Missing Resources_pyside.py: Run script update_resources.sh"
except ImportError:
    print "Can't import PySide modules."
    sys.exit()

import GetDistPlots
import MCSamples

import batchJob

# ==============================================================================

class MainWindow(QMainWindow):

    def __init__(self, app):
        """
        Initialize of GUI components.
        """
        super(MainWindow, self).__init__()

        self.app = app

        # GUI setup
        self.createWidgets()
        self.createActions()
        self.createMenus()
        self.createStatusBar()

        self.setAttribute(Qt.WA_DeleteOnClose)
        self.setWindowTitle("GetDist GUI")

        # Allow to shutdown the GUI with Ctrl+C
        signal.signal(signal.SIGINT, signal.SIG_DFL)

        # Path for .ini file
        self.iniFile = ""

        # Path of root directory
        self.rootdirname = None

        # Dict for roots
        self.rootnames = {}

        self._resetGridData()
        self._resetPlotData()


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

        self.helpMenu = self.menuBar().addMenu("&Help")
        self.helpMenu.addAction(self.aboutAct)

    def createStatusBar(self):
        """
        Create Qt status bar.
        """
        self.statusBar().showMessage("Ready", 2000)

    def createWidgets(self):
        """
        Create widgets.
        """
        self.centralWidget = QWidget(self)
        self.setCentralWidget(self.centralWidget)

        self.selectWidget = QWidget(self.centralWidget)
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
        self.listRoots.setDragEnabled(True)
        self.listRoots.viewport().setAcceptDrops(True)
        self.listRoots.setDropIndicatorShown(True)
        self.listRoots.setDragDropMode(QAbstractItemView.InternalMove)
        self.connect(self.listRoots,
                     SIGNAL("itemChanged(QListWidgetItem *)"),
                     self.updateListRoots)

        self.pushButtonRemove = QPushButton(QIcon(":/images/file_remove.png"),
                                            "", self.selectWidget)
        self.pushButtonRemove.setToolTip("Remove a root directory")
        self.connect(self.pushButtonRemove, SIGNAL("clicked()"),
                     self.removeRoot)

        self.comboBoxParamTag = QComboBox(self.selectWidget)
        self.comboBoxParamTag.clear()
        self.connect(self.comboBoxParamTag,
                     SIGNAL("activated(const QString&)"), self.setParamTag)

        self.comboBoxDataTag = QComboBox(self.selectWidget)
        self.comboBoxDataTag.clear()
        self.connect(self.comboBoxDataTag,
                     SIGNAL("activated(const QString&)"), self.setDataTag)

        self.comboBoxRootname = QComboBox(self.selectWidget)
        self.comboBoxRootname.clear()
        self.connect(self.comboBoxRootname,
                     SIGNAL("activated(const QString&)"), self.setRootname)

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
        layoutTop.addWidget(self.comboBoxRootname,  3, 0, 1, 4)
        layoutTop.addWidget(self.comboBoxParamTag,  3, 0, 1, 4)
        layoutTop.addWidget(self.comboBoxDataTag,   4, 0, 1, 4)
        layoutTop.addWidget(self.selectAllX,        5, 0, 1, 2)
        layoutTop.addWidget(self.selectAllY,        5, 2, 1, 2)
        layoutTop.addWidget(self.listParametersX,   6, 0, 5, 2)
        layoutTop.addWidget(self.listParametersY,   6, 2, 1, 2)
        layoutTop.addWidget(self.toggleFilled,      7, 2, 1, 2)
        layoutTop.addWidget(self.toggleLine,        8, 2, 1, 2)
        layoutTop.addWidget(self.toggleColor,       9, 2, 1, 1)
        layoutTop.addWidget(self.comboBoxColor,     9, 3, 1, 1)
        layoutTop.addWidget(self.trianglePlot,     10, 2, 1, 2)
        layoutTop.addWidget(self.pushButtonPlot,   12, 0, 1, 4)
        self.selectWidget.setLayout(layoutTop)

        self.listRoots.hide()
        self.pushButtonRemove.hide()
        self.comboBoxRootname.hide()
        self.comboBoxParamTag.hide()
        self.comboBoxDataTag.hide()

        self.plotWidget = QWidget(self.centralWidget)
        layout = QVBoxLayout(self.plotWidget)
        self.plotWidget.setLayout(layout)

        splitter = QSplitter(self.centralWidget)
        splitter.addWidget(self.selectWidget)
        splitter.addWidget(self.plotWidget)
        w = self.width()
        splitter.setSizes([w/4., 3*w/4.])

        hbox = QHBoxLayout()
        hbox.addWidget(splitter)
        self.centralWidget.setLayout(hbox)
        self.readSettings()

    def setIniFile(self, iniFile):
        logging.debug("ini file is %s"%iniFile)
        self.iniFile = iniFile

    def closeEvent(self, event):
        self.writeSettings()
        event.accept()

    def readSettings(self):
        settings = QSettings("cosmologist", "cosmomc_gui")
        pos = settings.value("pos", QPoint(200, 200))
        size = settings.value("size", QSize(400, 400))
        self.resize(size)
        self.move(pos)

    def writeSettings(self):
        settings = QSettings("cosmologist", "cosmomc_gui")
        settings.setValue("pos", self.pos())
        settings.setValue("size", self.size())

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
            self, "Choose a file name", '.', "Python (*.py)")
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
        rootname = None
        item = self.listRoots.currentItem()
        if item is not None:
            rootname = str(item.text())

        if rootname is None:
            QMessageBox.warning(self, "Marge Stats", "No rootname. Can't show marge stats")
            #logging.warning("No rootname. Can't show marge stats")
            return

        if self.plotter is None:
            QMessageBox.warning(self, "Marge Stats", "No plotter. Can't show marge stats")
            #logging.warning("No plotter. Can't show marge stats")
            self.statusBar().showMessage("No data available", 2000)
            return

        # FIXME: wait cursor
        text = self.plotter.sampleAnalyser.getMargeStats(rootname)
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
	settings = QSettings('cosmomc', 'gui')
        last_dir = settings.value('lastSearchDirectory')
        last_dir = ''
        if not last_dir: last_dir = os.getcwd()

        title = self.tr("Choose an existing case")
        dirName = QFileDialog.getExistingDirectory(
            self, title, last_dir,
            QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks)
        dirName = str(dirName)
        logging.debug("dirName: %s"%dirName)
        if dirName is None or dirName=='':
            return # No directory selected

        settings.setValue('lastSearchDirectory', dirName)

        # Check if it's a grid
        filebatch = os.path.join(dirName, "batch.pyobj")
        if os.path.isfile(filebatch):
            self.rootdirname = dirName
            self.lineEditDirectory.setText(self.rootdirname)
            self._readGridChain(self.rootdirname)
            return
        else:
            self._resetGridData()

        self.rootdirname = dirName
        self.lineEditDirectory.setText(self.rootdirname)

        # Update combo box
        filesparam = MCSamples.GetParamNamesFiles(self.rootdirname)
        self._updateComboBoxRootname(filesparam)

        if self.plotter is None:
            self.plotter = GetDistPlots.GetDistPlotter(
                mcsamples=True,
                ini_file=self.iniFile)

        # Hide combo boxes and fill list
        self.comboBoxParamTag.hide()
        self.comboBoxDataTag.hide()
        self.listRoots.show()
        self.pushButtonRemove.show()
        self.comboBoxRootname.show()

    def _updateParameters(self):
        all_params = []

        for k, values in self.rootnames.items():
            rootname, state = values
            if state:
                params = self.plotter.sampleAnalyser.usedParamsForRoot(rootname)
                logging.debug("%i parameters"%len(params))
                all_params.append(params)

        if all_params:
            if len(all_params)==1:
                paramNames = all_params[0]
                logging.debug("%i parameters used"%len(paramNames))
            else:
                # Find shared parameters
                shared_params = []
                for param in all_params[0]:
                    shared = True
                    for other_params in all_params[1:]:
                        if not param in other_params: shared = False
                    if shared: shared_params.append(param)
                paramNames = shared_params
                logging.debug("%i (shared) parameters used"%len(paramNames))

            self._updateListParametersX(paramNames)
            self._updateListParametersY(paramNames)
            self._updateComboBoxColor(paramNames)
        else:
            logging.warning("No parameters found")

    def _resetPlotData(self):
        # Plot parameters
        self.items_x = []
        self.items_y = []
        self.plotter = None

        # Script
        self.script = ""

    def _resetGridData(self):
        # Grid chains parameters
        self.is_grid = False
        self.grid_items = {}
        self.paramTag = ""
        self.dataTag = ""
        self.data2chains = {}

    def _readGridChain(self, batchPath):
        """
        Setup of a grid chain.
        """
        # Reset data
        self.rootnames = {}
        self._resetGridData()
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
        self.comboBoxRootname.hide()
        self.comboBoxParamTag.show()
        self.comboBoxDataTag.show()
        self.listRoots.show()
        self.pushButtonRemove.show()
        self._updateComboBoxParamTag(self.grid_items.keys())

    def _updateComboBoxRootname(self, listOfRoots):
        self.comboBoxParamTag.hide()
        self.comboBoxRootname.show()
        self.comboBoxRootname.clear()
        baseRoots = [ os.path.basename(root) for root in listOfRoots ]
        self.comboBoxRootname.addItems(baseRoots)

    def setRootname(self, strParamName):
        """
        Slot function called on change of comboBoxRootname.
        """
        paramName = str(strParamName)
        basename = paramName.replace('.paramnames', '')
        root = os.path.join(self.rootdirname, basename)

        if self.plotter is None:
            QMessageBox.warning(self, "Set root", "No plotter defined")
            return

        if not self.rootnames.has_key(basename):
            # FIXME: wait cursor
            #self.app.setOverrideCursor(QCursor(Qt.WaitCursor))
            self.plotter.sampleAnalyser.addRoot(root)
            logging.debug("Add root %s"%basename)
            self.rootnames[basename] = [root, True]
            item = QListWidgetItem(self.listRoots)
            item.setText(basename)
            item.setCheckState(Qt.Checked)
            idx = self.listRoots.count()
            self.listRoots.insertItem(idx, item)
            self._updateParameters()
            #self.app.restoreOverrideCursor()
        else:
            logging.warning("Root allready in list")

    def updateListRoots(self, item):
        logging.debug("updateListRoots")
        nitems = self.listRoots.count()
        for i in range(nitems):
            item = self.listRoots.item(i)
            root = item.text()
            state = item.checkState()==Qt.Checked
            if self.rootnames.has_key(root):
                self.rootnames[root][1] = state
            else:
                logging.warning('No root found for %s'%root)
        self._updateParameters()

    def removeRoot(self):
        logging.debug("Remove root")
        for i in range(self.listRoots.count()):
            item = self.listRoots.item(i)
            if item and item.isSelected():
                root = str(item.text())
                logging.debug("Remove root %s"%root)
                self.plotter.sampleAnalyser.removeOtherRoot(root)
                if self.rootnames.has_key(root):
                    del self.rootnames[root]
                self.listRoots.takeItem(i)
                self._updateParameters()

    def getRoots(self):
        logging.debug("Get status for roots")
        for i in range(self.listRoots.count()):
            item = self.listRoots.item(i)
            root = str(item.text())
            state = (item.checkState()==Qt.Checked)
            if self.rootnames.has_key(root):
                self.rootnames[root][1] = state
                logging.debug("status for '%s': %s"%(root, state))
            else:
                logging.warning("No key for root '%s'"%root)

    def _updateComboBoxParamTag(self, listOfParams=[]):
        if self.rootdirname and os.path.isdir(self.rootdirname):
            self.comboBoxRootname.hide()
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
        if self.rootnames.has_key(basename):
            logging.warning("Root name '%s' allready selected"%basename)
            return

        if self.plotter is None:
            self.plotter = GetDistPlots.GetDistPlotter(mcsamples=True,
                                                       ini_file=self.iniFile)
        self.statusBar().showMessage("Read chains in %s"%str(rootname))
        # FIXME: wait cursor
        self.plotter.sampleAnalyser.addRoot(rootname)
        self.rootnames[basename] = [rootname, True]
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
        items_x_old = self.items_x
        self.items_x = []
        self.listParametersX.clear()
        for item in items:
            listItem = QListWidgetItem()
            listItem.setText(item)
            listItem.setFlags(listItem.flags() | Qt.ItemIsUserCheckable)
            listItem.setCheckState(Qt.Unchecked)
            self.listParametersX.addItem(listItem)

        if items_x_old:
            for item_x in items_x_old:
                match_items = self.listParametersX.findItems(item_x, Qt.MatchExactly)
                if match_items:
                    logging.debug("Re check param %s"%item_x)
                    match_items[0].setCheckState(Qt.Checked)

    def itemParamXChanged(self, item):
        """
        Slot function called when item in listParametersX changes.
        """
        item_x = str(item.text())
        if item.checkState()==Qt.Checked:
            if item_x not in self.items_x:
                logging.debug("Add item %s"%item_x)
                self.items_x.append(item_x)
        else:
            if item_x in self.items_x:
                logging.debug("Del item %s"%item_x)
                self.items_x.remove(item_x)

    def _updateListParametersY(self, items):
        """
        """
        logging.debug("Fill ListParametersY with %i items"%(len(items)))
        items_y_old = self.items_y
        self.items_y = []
        self.listParametersY.clear()
        for item in items:
            listItem = QListWidgetItem()
            listItem.setText(item)
            listItem.setFlags(listItem.flags() | Qt.ItemIsUserCheckable)
            listItem.setCheckState(Qt.Unchecked)
            self.listParametersY.addItem(listItem)

        if items_y_old:
            for item_y in items_y_old:
                match_items = self.listParametersY.findItems(item_y, Qt.MatchExactly)
                if match_items:
                    logging.debug("Re check param %s"%item_y)
                    match_items[0].setCheckState(Qt.Checked)

    def itemParamYChanged(self, item):
        """
        Slot function called when item in listParametersY changes.
        """
        item_y = str(item.text())
        if item.checkState()==Qt.Checked:
            if item_y not in self.items_x:
                logging.debug("Add item %s"%item_y)
                self.items_y.append(item_y)
        else:
            if item_y in self.items_y:
                logging.debug("Del item %s"%item_y)
                self.items_y.remove(item_y)

    def statusSelectAllX(self):
        """
        Slot function called when selectAllX is modified.
        """
        self.items_x = []
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
        self.items_y = []
        if self.selectAllY.isChecked():
            state = Qt.Checked
        else:
            state = Qt.Unchecked
        for i in range(self.listParametersY.count()):
            self.listParametersY.item(i).setCheckState(state)

    def _updateComboBoxColor(self, listOfParams=[]):
        if self.rootdirname and os.path.isdir(self.rootdirname):
            param_old = str(self.comboBoxColor.currentText())
            self.comboBoxColor.clear()
            if not listOfParams:
                listOfParams = [
                    d for d in os.listdir(self.rootdirname)
                    if os.path.isdir(os.path.join(self.rootdirname, d)) ]
            self.comboBoxColor.addItems(listOfParams)

            idx = self.comboBoxColor.findText(param_old, Qt.MatchExactly)
            if idx<>-1:
                self.comboBoxColor.setCurrentIndex(idx)
                logging.debug("Re set param %s at index %i"%(param_old, idx))
            else:
                logging.debug("Param %s not found in new list"%(param_old))

    def plotData(self):
        """
        Slot function called when pushButtonPlot is pressed.
        """

        # Ensure at least 1 root name specified
        self.getRoots()

        nroots = 0
        for basename, values in self.rootnames.items():
            rootname, state = values
            if state: nroots += 1

        if nroots==0:
            logging.warning("No rootname selected")
            QMessageBox.warning(self, "Plot data", "No root selected")
            return

        if self.plotter is None:
            QMessageBox.warning(self, "Plot data", "No GetDistPlotter instance")
            return

        self.plotter.settings.setWithSubplotSize(3.000000)
        if self.plotter.fig is not None:
            self.plotter.fig.clf()

        plt.close('all')

        # X and Y items
        items_x = self.items_x
        items_y = self.items_y

        # Script
        self.script  = ""
        self.script += "import GetDistPlots, os\n"
        self.script += "g=GetDistPlots.GetDistPlotter(mcsamples=True, ini_file=r'%s')\n"%self.iniFile
        self.script += "g.settings.setWithSubplotSize(3.000000)\n"
        self.script += "outdir='.'\n"
        self.script += "dict_roots={}\n"

        roots = []

        # Root names
        for i in range(self.listRoots.count()):
            item = self.listRoots.item(i)
            basename = str(item.text())
            if item.checkState()==Qt.Checked:
                if self.rootnames.has_key(basename):
                    rootname, state = self.rootnames[basename]
                    roots.append(os.path.basename(rootname))
                    self.script += "dict_roots['%s'] = r'%s'\n"%(basename, rootname)
                else:
                    logging.warning("self.rootnames has no key %s"%basename)
            else:
                logging.debug("root '%s' not selected"%basename)

        self.script += "g.sampleAnalyser.addRoots(dict_roots.values())\n"
        self.script += "roots=dict_roots.keys()\n"

        logging.debug("Plotting with roots = %s"%str(roots))

        # Plot parameters
        filled = self.toggleFilled.isChecked()
        line = self.toggleLine.isChecked()
        color = self.toggleColor.isChecked()
        color_param = str(self.comboBoxColor.currentText())

        # Check type of plot
        if self.trianglePlot.isChecked():
            # Triangle plot
            if len(items_x)>1:
                params = items_x
                logging.debug("Triangle plot with params = %s"%str(params))
                self.script += "params = %s\n"%str(params)
                if color:
                    param_3d = color_param
                    self.script += "param_3d = '%s'\n"%str(color_param)
                else:
                    param_3d = None
                    self.script += "param_3d = None\n"
                self.script += "filled = %s\n"%filled
                try:
                    self.plotter.triangle_plot(roots, params, plot_3d_with_param=param_3d, filled_compare=filled)
                except:
                    QMessageBox.critical(
                        self, "Triangle plot",
                        "Error for command:\n\ntriangle_plot(roots, params, plot_3d_with_param=param_3d, filled_compare=filled)\n\nwith roots=%s\nparams=%s\nparam_3d=%s\nfilled=%s\n"%(str(roots), str(params), str(param_3d), str(filled)))
                    return
                self.updatePlot()
                self.script += "g.triangle_plot(roots, params, plot_3d_with_param=param_3d, filled_compare=filled)\n"
            else:
                QMessageBox.warning(self, "Triangle plot", "Select more than 1 X parameter")

        elif len(items_x)>0 and len(items_y)==0:
            # 1D plot
            params = items_x
            logging.debug("1D plot with params = %s"%str(params))
            self.script += "params=%s\n"%str(params)
            try:
                self.plotter.plots_1d(roots, params=params)
            except:
                QMessageBox.critical(
                        self, "Plot 1D",
                        "Error for command:\n\nplots_1d(roots, params=params)\n\nwith roots=%s\nparams=%s\n"%(str(roots), str(params)))
                return
            self.updatePlot()
            self.script += "g.plots_1d(roots, params=params)\n"

        elif len(items_x)>0 and len(items_y)>0:
            if filled or line:
                if len(items_x)>1 and len(items_y)>1:
                    # Rectangle plot
                    self.script += "xparams = %s\n"%str(items_x)
                    self.script += "yparams = %s\n"%str(items_y)
                    self.script += "filled=%s\n"%filled
                    logging.debug("Rectangle plot with xparams=%s and yparams=%s"%(str(items_x), str(items_y)))
                    try:
                        self.plotter.rectangle_plot(items_x, items_y, roots=roots)
                    except:
                        QMessageBox.critical(
                            self, "Plot 2D",
                            "Error for command:\n\nrectangle_plot(xparams, yparams, roots=roots)\n\nwith xparams=%s\nyparams=%s\nroots=%s\n"%(str(items_x), str(items_y), str(roots)))
                        return
                    self.updatePlot()
                    self.script += "g.rectangle_plot(xparams, yparams, roots=roots)\n"

                else:
                    # 2D plot
                    if len(items_x)==1 and len(items_y)==1:
                        pairs = [ [items_x[0], items_y[0]] ]
                    elif len(items_x)==1 and len(items_y)>1:
                        item_x = items_x[0]
                        pairs = zip([item_x] * len(items_y), items_y)
                    elif len(items_x)>1 and len(items_y)==1:
                        item_y = items_y[0]
                        pairs = zip(items_x, [item_y] * len(items_x))
                    else:
                        pairs = []
                    self.script += "pairs = %s\n"%pairs
                    logging.debug("2D plot with pairs = %s"%str(pairs))
                    self.script += "filled=%s\n"%filled
                    try:
                        self.plotter.plots_2d(roots, param_pairs=pairs, filled=filled)
                    except:
                        QMessageBox.critical(
                            self, "Plot 2D",
                            "Error for command:\n\nplots_2d(roots, param_pairs=pairs, filled=filled)\n\nwith roots=%s\npairs=%s\nfilled=%s\n"%(str(roots), str(pairs), str(filled)))
                        return
                    self.updatePlot()
                    self.script += "g.plots_2d(roots, param_pairs=pairs, filled=filled)\n"

            if color:
                # 3D plot
                sets = []
                sets.append(items_x)
                x = items_x[0]
                y = items_y[0]
                sets = [[x, y, color_param]]
                logging.debug("3D plot with sets = %s"%str(sets))
                self.script += "sets = []\n"
                self.script += "sets.append(['%s', '%s', '%s'])\n"%(x, y, color_param)
                try:
                    self.plotter.plots_3d(roots, sets)
                except:
                    QMessageBox.critical(
                        self, "Plot 3D",
                        "Error for command:\n\nplots_3d(roots, sets)\n\nwith roots=%s\nsets=%s\n"%(str(roots), str(sets)))
                    return
                self.updatePlot()
                self.script += "g.plots_3d(roots, sets)\n"

        else:
            text  = ""
            text += "Wrong parameters selection. Specify parameters such as:\n"
            text += "\n"
            text += "Triangle plot: Click on 'Triangle plot' and select more than 1 X parameters\n"
            text += "\n"
            text += "1D plot: Select X parameter(s)\n"
            text += "\n"
            text += "2D plot: Select X parameter(s), Y parameter(s) and select 'Filled' or 'Line'\n"
            text += "\n"
            text += "3D plot: Select X parameter, Y parameter and 'Color by' parameter\n"
            text += "\n"

            QMessageBox.warning(self, "Plot usage", text)
            self.script = ""
            return

        self.script += "g.export(os.path.join(outdir, 'plot.pdf'))\n"


    def updatePlot(self):
        if self.plotter.fig is None:
            #logging.debug("Define an empty central widget")
            #self.figure = None
            self.canvas = None
        else:
            #logging.debug("Define new canvas in central widget")
            # Remove everything from layout
            i = 0
            while 1:
                item = self.plotWidget.layout().takeAt(i)
                if item is None: break
            if hasattr(self, "canvas"): del self.canvas
            if hasattr(self, "toolbar"): del self.toolbar

            self.canvas = FigureCanvas(self.plotter.fig)
            self.toolbar = NavigationToolbar(self.canvas, self)
            self.plotWidget.layout().addWidget(self.toolbar)
            self.plotWidget.layout().addWidget(self.canvas)
            self.canvas.draw()
            #self.plotWidget.update()
            self.plotWidget.show()

# ==============================================================================

class DialogMargeStats(QDialog):

    def __init__(self, parent=None, text=""):
        QDialog.__init__(self, parent)

        self.label = QLabel(self)
        self.table = QTableWidget(self)
        self.table.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.table.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        #self.table.horizontalHeader().hide()
        self.table.verticalHeader().hide()

        layout = QGridLayout()
        layout.addWidget(self.label, 0, 0)
        layout.addWidget(self.table, 1, 0)
        self.setLayout(layout)

        self.setWindowTitle(self.tr("Dialog for MargeStats"))

        if (text):
            lines = text.split("\n")
            line0 = lines.pop(0)
            self.label.setText(line0)
            line = lines.pop(0) # empty
            line = lines.pop(0) # headers
            headers = line.split(' ')
            headers = [ h for h in headers if h<>'' ] + [ 'Label' ]
            self.table.setColumnCount(len(headers))
            self.table.setHorizontalHeaderLabels(headers)
            self.table.verticalHeader().setVisible(False)
            self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
            self.table.setSelectionMode(QAbstractItemView.SingleSelection)
            self.table.setRowCount(len(lines))
            irow = 0
            for line in lines:
                values = line.split(' ')
                values = [ v for v in values if v<>'' ]
                icol = 0
                for value in values:
                    item = QTableWidgetItem(value)
                    item.setFlags(item.flags() ^ Qt.ItemIsEditable)
                    self.table.setItem(irow, icol, item)
                    icol += 1
                irow += 1

            self.table.resizeRowsToContents()
            self.table.resizeColumnsToContents()

            w = self.table.horizontalHeader().length() + 40
            h = self.table.verticalHeader().length()   + 40
            self.resize(w, h)

# ==============================================================================

if __name__ == "__main__":

    app = QApplication(sys.argv)
    mainWin = MainWindow(app)
    mainWin.show()
    sys.exit(app.exec_())

# ==============================================================================


