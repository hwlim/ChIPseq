#!/usr/bin/env python

import sys, os, string, glob
from os import path
import pipeUtil
import genomeUtil
from multiprocessing import Pool
import BIN2NRM


def Job((srcDir, desDir, binsize)):
##    print srcDir
##    print desDir
##    print binsize
    totalTagCnt = BIN2NRM.BIN2NRM_DIR(srcDir, desDir, binsize)


## sys.argv[1]: Max number of processes
## sys.argv[2]: desDir, (default: NRM)
## sys.argv[3]: srcDir, (default: BIN)
if __name__=='__main__':
    if len(sys.argv)<4:
        print('Usage:\n\tS3.BIN2NRM [MaxProc] [desDir] [srcDir]')
        sys.exit(1)

    srcDir = sys.argv[3]
    desDir = sys.argv[2]
    MaxProc = int(sys.argv[1])

    
    genome, binsize, clonal_mode, fragLen, shiftMode, fgD, namePairL, fragLenD = \
            pipeUtil.readConfig(LEVEL=1)
##    chrLenD = genomeUtil.getChrLenD(genome)

    os.system('mkdir -p %s' % desDir)

    fgL = reduce(lambda x,y: x+y, fgD.values())
    print fgL
##    print fgL
##    binDirL = [path.join(srcDir, fg) for fg in fgL] #glob.glob(path.join(tsvDir, '*'))
##    nrmDirL = [binDir.replace(srcDir, desDir) for binDir in binDirL]
##    argL = [[path.join(tsvDir, fg), path.join(binDir, fg), \
##             clonal_mode, binsize, chrLenD, fragLenD[fg]] \
##            for fg in fgL]

    
##    binDirL = glob.glob(path.join(srcDir, '*'))
##    desDirL = [binDir.replace(srcDir, desDir) for binDir in binDirL]
    argL = [[path.join(srcDir,fg), path.join(desDir,fg), binsize] \
            for fg in fgL]
    
    pool = Pool(processes=MaxProc)
    try:
        pool.map(Job, argL)
    except KeyboardInterrupt:
        pool.terminate()

        

##
##
##genome, binsize, MaxProc, monoclonal, fgL, bgL = pipeUtil.readConfig()
##if len(bgL)==0: bgD = None
##else: bgD = dict([fgL[i], bgL[i]] for i in xrange(len(fgL)))
##        
##homerDir = "HomerTmp/"
##nrmDir = "nrm/"
##os.system('mkdir -p nrm')
##
##All = []
##for dataName in fgL+bgL:
##        All = All + Util.GetDirList(homerDir+dataName+"/*.bin",DIR=1)
##        os.system('mkdir -p %s%s' % (nrmDir, dataName))
##
##for binFile in All:
##        dataName = binFile.split("/")[-2]
##        nrmName = binFile.replace(homerDir, nrmDir).replace(".bin",".nrm")
##        divName = nrmName.replace(".nrm", ".div")
##
##        H = Util.ReadColumn(binFile,[0])
##
##        InfoFile = os.path.dirname(binFile)+"/tagInfo.txt"
##        total, fragment = GetInfo(InfoFile)
##
##        NormFG = RPKM(H,total)
##        PrintHistone(NormFG,nrmName)
##
##
##			
