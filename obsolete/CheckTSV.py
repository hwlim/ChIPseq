#!/usr/bin/env python

import os, sys, string, glob
from os import path
import pipeUtil

if len(sys.argv)<2: srcDir = './TSV'
else: srcDir = sys.argv[1]

hDirL = [d for d in glob.glob(path.join(srcDir,'*')) if os.path.isdir(d)]
##hDirL = os.listdir(srcDir)
##print hDirL

fout = open('checkTSV.info', 'w')
fout.write('DataName\tTotalTagCnt\tFragLength\tAvgTagPerPos\tUniqPositions\n')
for hDir in sorted(hDirL):
    dataName=hDir.split('/')[-1]
    TTC, fLen, avgTPP, Position = pipeUtil.readTSVInfo(hDir)
    fout.write('%s\t%d\t%d\t%.3f\t%d\n' %(dataName, TTC, fLen, avgTPP, Position))
fout.close()
