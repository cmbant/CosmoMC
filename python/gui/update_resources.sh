#!/bin/bash

# If Qt module is PyQt4
pyrcc4 -o Resources_pyqt.py Resources.qrc 

# If Qt module is PySide
pyside-rcc -o Resources_pyside.py Resources.qrc 

