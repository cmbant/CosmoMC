# -*- coding: utf-8 -*-

import os
import sys
import logging

try: import argparse
except:
    print 'use "module load" to load python 2.7'
    sys.exit()

try:
    from PySide.QtGui  import QApplication
    os.environ['QT_API'] = 'pyside'
except ImportError:
    print "Can't import PySide modules."
    sys.exit()

from MainWindow import MainWindow

parser = argparse.ArgumentParser(description='GetDist GUI')
parser.add_argument('-v', '--verbose', help='verbose', action="store_true")
parser.add_argument('--ini', help='Path to .ini file', default=None)
args = parser.parse_args()

# Configure the logging
level = logging.INFO
if args.verbose: level = logging.DEBUG
FORMAT = '%(asctime).19s [%(levelname)s]\t[%(filename)s:%(lineno)d]\t\t%(message)s'
logging.basicConfig(level=level, format=FORMAT)

# Change to suitable directory
cosmomc_dir = os.path.normpath(os.path.dirname(os.path.abspath(__file__)) + '/../../') + os.sep

# GUI application
app = QApplication(sys.argv)
mainWin = MainWindow(app, cosmomc_dir, args.ini)
mainWin.show()
sys.exit(app.exec_())
