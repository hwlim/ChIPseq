#!/usr/bin/env python
#$ -S /gpfs/fs0/share/python2.6.5/Python-2.6.5/python

import sys, string, os, glob
from os import path
import pipeUtil
##import genomeUtil
import BIN2WIG


if __name__ == '__main__':
    if len(sys.argv)<3:
        print 'Arguments: desDir, srcDir, (viewLimits), (height)'
        sys.exit(1)

    genome, binsize, clonal_mode, fragLen, shiftMode, fgD, namePairL, fragLenD = \
            pipeUtil.readConfig(LEVEL=1)

    desDir, srcDir = \
            path.normpath(sys.argv[1]), path.normpath(sys.argv[2])
    
    if len(sys.argv)<5: height=None
    else: height = int(sys.argv[4])
    if len(sys.argv)<4: viewLimits=None
    else: viewLimits = int(sys.argv[3])
    os.system('mkdir -p %s' % desDir)
    
    for (src, name) in sorted(namePairL):
        if not src in fgD:
            print 'ToWig %s --> %s' % (path.join(srcDir,src), desDir+'/'+src+'.wig')
            BIN2WIG.BIN2WIG(path.join(srcDir,src), desDir, binsize, \
                            name, viewLimits, height)
            dataName = path.basename(srcDir)
            os.system('gzip %s' % path.join(desDir, src+'.wig'))
            
