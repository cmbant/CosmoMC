import sys,re

with open(sys.argv[1],'r') as f:
    f.readline()
    f.readline()
    f.readline()
    x=f.readline().strip()
    while(not x==''):
#            print(x.split())
            print("param["+x.split()[2]+"]="+x.split()[1])
            x=f.readline().strip()



