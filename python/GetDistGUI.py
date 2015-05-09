# -*- coding: utf-8 -*-

import os
import sys
import logging

try:
    import argparse
except ImportError:
    print 'use "module load" to load python 2.7, or see docs/readme_python.html for how to install'
    sys.exit()

try:
    from PySide.QtGui import QApplication
    os.environ['QT_API'] = 'pyside'
except ImportError:
    print "Can't import PySide modules. See docs/readme_python.html for how to install."
    sys.exit()

try:
    from getdist.gui.mainwindow import MainWindow
except ImportError:
    print "Configure your PYTHONPATH as described in the readme!"
    sys.exit()

parser = argparse.ArgumentParser(description='GetDist GUI')
parser.add_argument('-v', '--verbose', help='verbose', action="store_true")
parser.add_argument('--ini', help='Path to .ini file', default=None)
args = parser.parse_args()

# Configure the logging
level = logging.INFO
if args.verbose:
    level = logging.DEBUG
FORMAT = '%(asctime).19s [%(levelname)s]\t[%(filename)s:%(lineno)d]\t\t%(message)s'
logging.basicConfig(level=level, txformat=FORMAT)

# GUI application
app = QApplication(sys.argv)
mainWin = MainWindow(app, ini=args.ini)
mainWin.show()
sys.exit(app.exec_())
