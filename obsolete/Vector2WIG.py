#!/usr/bin/env python
# temporary, auxilary utility
import glob, os, sys, string
from os import path
import util
import genomeUtil


# Arguments
# src, des, Chr, start, step, viewLimits, name, description
if __name__=='__main__':
    NofArg = len(sys.argv)
    if len(sys.argv) < 4:
        print 'Arguments: src, des, Chr, start, step, viewLimits, name, description'
        sys.exit(1)

    srcFile = sys.argv[1]
    desFile = sys.argv[2]
    chrStr = sys.argv[3]

    if NofArg < 5: start = 1
    else: start=int(sys.argv[4])
    if NofArg < 6: step = 100
    else: step=int(sys.argv[5])
    if NofArg < 7: viewLimits = None
    else: viewLimits=int(sys.argv[6])
    if NofArg < 8: name = None
    else: name=sys.argv[7]
    if NofArg < 9: description = name
    else: description=sys.argv[8]
    if viewLimits==0: viewLimits==None
    
    if not path.exists(srcFile):
        raise RuntimeError("Source file does not exist: %s" % srcFile)

##    data = util.readOneCol(srcFile, str)
    genomeUtil.writeWigHeader(desFile, chrStr, start, step, viewLimits, name, description)
    os.system('cat %s >> %s' % (srcFile, desFile))
