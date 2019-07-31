#!/usr/bin/env python

import sys, string, os, glob
from os import path
import pipeUtil
import Util
from multiprocessing import Pool

##srcDir = "BED"
##desDir = "TSV"


def BED2TSV((srcFile, desDir, genome)):
    dataName = path.basename(srcFile).replace('.bed', '')
##    print dataName
    tsvDir = path.join(desDir, dataName)
    cmd = 'makeTagDirectory %s %s -format bed -genome %s -checkGC -normGC default' % (tsvDir, srcFile, genome)
    print cmd
    os.system(cmd)
    cmd = 'drawGCplot.sh %s' % (tsvDir)
    print cmd
    os.system(cmd)




# Command line argument: Number of process
if __name__ == '__main__':
##    bedL = glob.glob(path.join(srcDir, '*.bed'))
    if len(sys.argv)<4:
        print 'Usage:\n\tS1.BED2TSV_GC.py [MaxProc] [desDir] [srcDir]'
        sys.exit(1)
##        MaxProc = int(raw_input("Enter the number of processes: "))
    MaxProc = int(sys.argv[1])
    desDir = sys.argv[2]
    srcDir = sys.argv[3]

    genome, binsize, clonal_mode, fragLen, shiftMode, fgD, nameD, fragLenD = \
            pipeUtil.readConfig(LEVEL=1)
    
    os.system('mkdir -p %s' % desDir)

    fgL = reduce(lambda x,y: x+y, fgD.values())
    bedL = [path.join(srcDir, fg)+'.bed' for fg in fgL]
##    bedL = [path.join(srcDir, fg)+'.bed' for fg in fgL] #glob.glob(path.join(tsvDir, '*'))
####    desDirL = [srcDir.replace(tsvDir, binDir) for srcDir in srcDirL]
##    argL = [[path.join(tsvDir, fg), path.join(binDir, fg), \
##             clonal_mode, binsize, chrLenD, fragLenD[fg]] \
##            for fg in fgL]
    
    argL = [[bedFile, desDir, genome] for bedFile in bedL]

    print 'Starting multiprocessing: %d' % MaxProc
    pool = Pool(processes=MaxProc)
    try:
        pool.map(BED2TSV, argL)
    except KeyboardInterrupt:
        pool.terminate()

