#!/usr/bin/env python
""" unzip.py
    Copyright 2011 Hin-Tak Leung
    All rights reserved.
"""

import sys
import zipfile
import os
import os.path
import getopt
from datetime import datetime
import time

class unzip:
    def __init__(self, info_only=False):
        self.info_only = info_only

    def extract(self, file):
        zf = zipfile.ZipFile(file, allowZip64=True)
        if (self.info_only):
            # zf.printdir() #debug code
            for info in zf.infolist():
                try:
                    if info.filename.endswith('/'):  # test for info.CRC?
                        # need to be 1!
                        print "Path =", info.filename.rsplit("/", 1)[0].decode("GB2312").encode("UTF8")
                        print "Folder = +"
                    else:
                        print "Path =", info.filename.decode("GB2312").encode("UTF8")
                        print "Folder = -"
                except UnicodeDecodeError:
                    print
                    print "Failed with GB2312 =", info.filename.decode("GB2312", 'ignore').encode("UTF8", 'ignore')
                    print "Retry with GB18030 =", info.filename.decode("GB18030", 'ignore').encode("UTF8", 'ignore')
                    pass
                print "Size =", info.file_size
                print "Packed Size =", info.compress_size
                print "Modified = %04d-%02d-%02d %02d:%02d:%02d" % info.date_time
                print info.compress_type
                print info.comment
                print info.extra
                print info.create_system
                print info.create_version
                print info.extract_version
                print info.reserved
                print info.flag_bits
                print info.volume
                print info.internal_attr
                print info.external_attr
                print info.header_offset
                print "CRC = %06X" % info.CRC
            return

        # extract() is new to python 2.6
        # extract files to directory structure
        dirs = []
        for i, info in enumerate(zf.infolist()):
            try:
                newname = info.filename.decode("GB2312").encode("UTF8")
            except:
                newname = info.filename.decode("GB18030").encode("UTF8")
                print "Failed with GB2312 =", info.filename.decode("GB2312", 'ignore').encode("UTF8", 'ignore')
                print "Retry with GB18030 =", info.filename.decode("GB18030", 'ignore').encode("UTF8", 'ignore')
                pass
            zf.extract(info)
            os.rename(info.filename, newname)
            timestamp = datetime(info.date_time[0], info.date_time[1], info.date_time[2],
                                 info.date_time[3], info.date_time[4], info.date_time[5])
            epoch = time.mktime(timestamp.timetuple())
            os.utime(newname, (epoch, epoch))
            if info.filename.endswith('/'):
                dirs.append(info.filename)
        # while is dir deletion needed?
        dirs.sort(reverse=True)
        for dir in dirs:
            try:
                os.rmdir(dir)
            except:
                print "rmdir failed:", dir.decode("GB18030", 'ignore').encode("UTF8", 'ignore')
                pass

def usage():
    print """usage: unzip.py -z <zipfile>
    <zipfile> is the source zipfile to extract

    -z zipfile to extract
    -i  - show content info only

    long options also work:
    --zipfile=<zipfile>"""


def main():
    shortargs = 'ihz:'
    longargs = ['info', 'help', 'zipfile=']

    unzipper = unzip()

    try:
        opts, args = getopt.getopt(sys.argv[1:], shortargs, longargs)
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    zipsource = ""

    for o, a in opts:
        if o in ("-z", "--zipfile"):
            zipsource = a
        if o in ("-i", "--info"):
            unzipper.info_only = True
        if o in ("-h", "--help"):
            usage()
            sys.exit()

    if zipsource == "":
        usage()
        sys.exit()

    unzipper.extract(zipsource)

if __name__ == '__main__': main()
