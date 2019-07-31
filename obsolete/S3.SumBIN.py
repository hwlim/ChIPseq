#!/usr/bin/env python
# Optional step: Sum BIN files for multiple replicates
# Need to change for a flat folder structure
#   not a Cell/Histone structure
import sys, string, os, glob
from os import path
import pipeUtil
import SumBIN

srcDir = 'BIN'
desDir = 'sBIN'

##if len(sys.argv)<2: dataMode = raw_input("Enter data set mode (ENCODE, EPIATLAS): ")
##else: dataMode = sys.argv[1]

os.system('mkdir -p %s' % desDir)

binDirL = glob.glob(path.join(srcDir, '*'))

dirD = {}
for binDir in binDirL:
    dirName = path.basename(binDir)
    c, h, serial = pipeUtil.splitName(dirName, dataMode)
    h = pipeUtil.normHistoneName(h) # meaningful only for ENCODE data
    if not dirD.has_key(c): dirD[c] = {}
    if not dirD[c].has_key(h): dirD[c][h] = []
    dirD[c][h].append(binDir)

for c in dirD:
    for h in dirD[c]:
        print 'SumBIN: %s, %s' % (c, h)
        print 'Source List: %s' % '\t'.join(dirD[c][h])
        SumBIN.SumBIN(dirD[c][h], path.join(desDir, c+'.'+h))
