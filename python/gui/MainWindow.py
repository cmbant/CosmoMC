# -*- coding: utf-8 -*-

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

from CentralWidget import CentralWidget

class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
 
        self.setAttribute(Qt.WA_DeleteOnClose)
 
        self.textEdit = QTextEdit()
        self.setCentralWidget(self.textEdit)
 
        self.createActions()
        self.createMenus()
        self.createToolBars()
        self.createStatusBar()

        cw = CentralWidget()
        self.setCentralWidget(cw)

        self.setWindowTitle("GetDist GUI")
        self.resize(400, 300)

    def save(self):
        filename, filtr = QFileDialog.getSaveFileName(
            self,
            "Choose a file name", '.', "PNG (*.png)")
        if not filename:
            return

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

    def createActions(self):
        self.saveAct = QAction(
            "&Save", self,
            shortcut=QKeySequence.Save,
            statusTip="Save image as PNG", triggered=self.save)

        self.exportAct = QAction(
            "&Export", self,
            statusTip="Export image as PDF", triggered=self.export)

        self.scriptAct = QAction(
            "Script", self,
            statusTip="Export commands to script", triggered=self.script)

        self.exitAct = QAction(
            "E&xit", self, shortcut="Ctrl+Q",
            statusTip="Exit application",
            triggered=qApp.closeAllWindows)
 
        self.aboutAct = QAction(
            "&About", self,
            statusTip="Show About box",
            triggered=self.about)

    def createMenus(self):
        self.fileMenu = self.menuBar().addMenu("&File")
        self.fileMenu.addAction(self.saveAct)
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
        self.fileToolBar.addAction(self.saveAct)
        self.fileToolBar.addAction(self.exportAct)
        self.fileToolBar.addAction(self.exitAct)
 
    def createStatusBar(self):
        self.statusBar().showMessage("Ready", 2000)

