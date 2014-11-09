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
parser.add_argument('--ini', help='Path to .ini file', default='batch1/getdist_common.ini')
args = parser.parse_args()

# Configure the logging
level = logging.INFO
if args.verbose: level = logging.DEBUG
FORMAT = '%(asctime).19s [%(levelname)s]\t[%(filename)s:%(lineno)d]\t\t%(message)s'
logging.basicConfig(level=level, format=FORMAT)

# Change to suitable directory
os.chdir(os.path.dirname(os.path.abspath(__file__)) + '/../../')

# Ini file
iniFile = args.ini
if not os.path.isfile(iniFile):
    logging.warning("Invalid file %s"%iniFile)
    iniFile=""

# GUI application
app = QApplication(sys.argv)
mainWin = MainWindow(app)
if iniFile: mainWin.setIniFile(iniFile)
mainWin.show()
sys.exit(app.exec_())
