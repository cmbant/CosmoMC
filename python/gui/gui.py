# -*- coding: utf-8 -*-

import os
import sys

try: import argparse
except:
    print 'use "module load" to load python 2.7'
    sys.exit()

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


parser = argparse.ArgumentParser(description='GetDist GUI')

parser.add_argument('--ini', help='Path to .ini file', default='../batch1/getdist_common.ini')
args = parser.parse_args()

iniFile = args.ini
if not os.path.isfile(iniFile):
    print "Invalid file %s"%iniFile
    iniFile=""

app = QApplication(sys.argv)
mainWin = MainWindow()
if iniFile: 
    mainWin.setIniFile(iniFile)
mainWin.show()
sys.exit(app.exec_())
