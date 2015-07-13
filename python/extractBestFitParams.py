from __future__ import absolute_import
from __future__ import print_function
import sys

with open(sys.argv[1], 'r') as f:
    f.readline()
    f.readline()
    f.readline()
    x = f.readline().strip()
    while not x == '':
            print(("param[" + x.split()[2] + "]=" + x.split()[1]))
            x = f.readline().strip()



