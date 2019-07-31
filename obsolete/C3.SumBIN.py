#!/usr/bin/env python

import sys, string, os, glob
from os import path
import pipeUtil

srcDir = 'BIN'
desDir = 'sBIN'

if len(sys.argv)<2: dataMode = raw_input("Enter data set mode (ENCODE, EPIATLAS): ")
else: dataMode = sys.argv[1]

os.system('mkdir -p %s' % desDir)
os.system('rm ./log/*.*')

binDirL = glob.glob(path.join(srcDir, '*'))

dirD = {}
for binDir in binDirL:
    dirName = path.basename(binDir)
    c, h, serial = pipeUtil.splitName(dirName, dataMode)
    h = pipeUtil.normHistoneName(h)
    if not dirD.has_key(c): dirD[c] = {}
    if not dirD[c].has_key(h): dirD[c][h] = []
    dirD[c][h].append(binDir)

for c in dirD:
    for h in dirD[c]:
        jobName = c+'.'+h
        print 'Submitting a job: %s, %s' % (c, h)
        cmd = 'qsub -cwd -V -o log/%s.out -e log/%s.err ../scripts/SumBIN.py %s %s' \
              % (jobName, jobName, path.join(desDir, c, h), ' '.join(dirD[c][h]))
##        print cmd
        os.system(cmd)
