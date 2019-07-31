#!/usr/bin/env python

import sys, os, string, glob
from os import path
import pipeUtil
import genomeUtil
from multiprocessing import Pool
import TSV2BIN


##tsvDir = 'TSV'
##binDir = 'BIN'

##def Job((srcDir, desDir, clonal_mode, binsize, chrLenD, fragLen, shiftMode)):
##    print 1


def Job((srcDir, desDir, clonal_mode, binsize, chrLenD, fragLen, shiftMode)):
##    print srcDir
##    print desDir
##    print clonal_mode, binsize
##    os.system('mkdir -p %s' % desDir)
    totalTagCnt = TSV2BIN.TSV2BIN_DIR(srcDir, desDir, clonal_mode, binsize, chrLenD, fragLen, shiftMode)

##    fout = open(path.join(desDir, 'tagInfo.txt'), 'w')
##    fout.write('TotalTagCount=%d\n' % totalTagCnt)
##    fout.close()


## sys.argv[1]: Max number of processes
## sys.argv[2]: desDir, (default: BIN)
## sys.argv[3]: srcDir, (default: TSV)
if __name__=='__main__':
##    if len(sys.argv)<2: MaxProc = int(raw_input("Enter the number of processes: "))
##    else: MaxProc = int(sys.argv[1])

    if len(sys.argv)<4:
        print('Usage:\n\tS2.TSV2BIN [MaxProc] [desDir] [srcDir]')
        sys.exit(1)

    srcDir = sys.argv[3]
    desDir = sys.argv[2]
    MaxProc = int(sys.argv[1])

    

    genome, binsize, clonal_mode, fragLen, shiftMode, fgD, namePairL, fragLenD = \
            pipeUtil.readConfig(LEVEL=1)
    print 'Start converting TSV --> BIN'
    print 'Genome: %s' % genome
    print 'Bin Size: %d' % binsize
    print 'Clonal_mode: %d' % clonal_mode
    print 'Fragment length: %d' % fragLen
    print 'Shift mode: %d' % shiftMode
##    print fragLenD.items()
##    print fgD
    
    if fragLen==0: fragLen=None
    chrLenD = genomeUtil.getChrLenD(genome)

##    print fgD
    fgL = reduce(lambda x,y: x+y, fgD.values())
    print fgL
##    print fgL
##    srcDirL = [path.join(tsvDir, fg) for fg in fgL] #glob.glob(path.join(tsvDir, '*'))
##    desDirL = [srcDir.replace(tsvDir, binDir) for srcDir in srcDirL]
    argL = [[path.join(srcDir, fg), path.join(desDir, fg), \
             clonal_mode, binsize, chrLenD, fragLenD[fg], shiftMode] \
            for fg in fgL]
##    print argL
##    print MaxProc
##    argL = [[srcDirL[i], desDirL[i], clonal_mode, binsize, chrLenD, fragLen] \
##            for i in xrange(len(srcDirL))]
    
    pool = Pool(processes=MaxProc)
    try:
        pool.map(Job, argL)
##        pool.map(Job, range(1,5))
    except KeyboardInterrupt:
        pool.terminate()


