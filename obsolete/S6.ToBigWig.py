#!/usr/bin/env python
#$ -S /gpfs/fs0/share/python2.6.5/Python-2.6.5/python

import sys, string, os, glob
from os import path
##from optparse import OptionParser
import pipeUtil
##import genomeUtil


if __name__ == '__main__':
    if len(sys.argv)<3:
        print 'Usage: BIN2BW.sh [srcDir] [desFile] [binsize]'
        sys.exit(1)

    genome, binsize, clonal_mode, fragLen, shiftMode, fgD, namePairL, fragLenD = \
            pipeUtil.readConfig(LEVEL=1)

    desDir, srcDir = \
            path.normpath(sys.argv[1]), path.normpath(sys.argv[2])
    
    os.system('mkdir -p %s' % desDir)
    
    for (src, name) in sorted(namePairL):
        if not src in fgD:
            print 'ToBigWig %s --> %s' % (path.join(srcDir,src), path.join(desDir,name+'.bw'))
            FILEL=glob.glob(path.join(srcDir,src,'chr*.*'))
