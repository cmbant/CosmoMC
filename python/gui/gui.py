# -*- coding: utf-8 -*-

import sys

from MainWindow import MainWindow

try:
    from PyQt4.QtCore import *
    from PyQt4.QtGui  import *
except ImportError:
    try:
        from PySide.QtCore import *
        from PySide.QtGui  import *
    except ImportError:
        print "Can't import PyQt4.QtCore or PyQt4.QtGui modules." 
        sys.exit()

app = QApplication(sys.argv)
mainWin = MainWindow()
mainWin.show()
sys.exit(app.exec_())
