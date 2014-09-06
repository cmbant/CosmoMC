# -*- coding: utf-8 -*-

import os
import sys

try:
    from PyQt4.QtGui  import QApplication
    os.environ['QT_API'] = 'pyqt'
except ImportError:
    try:
        from PySide.QtGui  import QApplication
        os.environ['QT_API'] = 'pyside'
    except ImportError:
        print "Can't import PyQt4 or PySide modules." 
        sys.exit()

from MainWindow import MainWindow

app = QApplication(sys.argv)
mainWin = MainWindow()
mainWin.show()
sys.exit(app.exec_())
