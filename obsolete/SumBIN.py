#!/usr/bin/env python
#$ -S /gpfs/fs0/share/python2.6.5/Python-2.6.5/python
import sys, string, os, glob
from os import path
import pipeUtil
import genomeUtil


def SumBIN(dataDirL, desDir):    
    if len(dataDirL)<1: raise('At least one source data directory must be given.')

##    print dataDirL
    os.system('mkdir -p %s' % desDir)
    NF = len(dataDirL)
    genome = pipeUtil.readConfig(0)[0]
    NofChr = genomeUtil.getChrCnt(genome)
    dataFileL = [genomeUtil.getChrStr(chN, genome)+'.bin' \
                 for chN in xrange(1,NofChr+1)]  # except chrY

    for dataFile in dataFileL:
##        print 'Summing %s' % dataFile
        srcFileL = [path.join(dataDir, dataFile) for dataDir in dataDirL]
        gawkCmd = "'{if(NF!=%d){print %s; exit 1} sum=0; for(i=1;i<=NF;i++) sum+=$i; printf \"%%.3f\\n\", sum}'"\
                  % (NF, '"Error: Number of rows does not match"')
        cmd = 'paste %s | gawk %s > %s' \
              % (' '.join(srcFileL), gawkCmd, path.join(desDir, dataFile))
##        print cmd
        os.system(cmd)

    totalTagCnt = 0
    for dataDir in dataDirL:
##        print pipeUtil.readBINInfo(dataDir)
        totalTagCnt = totalTagCnt + pipeUtil.readBINInfo(dataDir)

    with open(path.join(desDir, 'tagInfo.txt'), 'w') as fout:
        fout.write('TotalTagCnt=%d\n' % totalTagCnt)
        fout.write('NumberOfRep=%d\n' % len(dataDirL))
        
        

# sys.argv[1]: destination Dir
# sys.argv[2:]: a list of source directory
# Example) SumBIN BIN/desDir BIN/srcDir1 BIN/srcDir2 ...
if __name__=='__main__':
    N = len(sys.argv)-2     # number of source data directory
    if N<1:
        print 'Usage:\n\tSumBIN.py [desDir] [a list of srcDir]'
        sys.exit(1)

    desDir = sys.argv[1]
    dataDirL = sys.argv[2:]
    for dataDir in dataDirL:
        if not path.isdir(dataDir):
            raise RuntimeError('Source directory does not exists: %s' % dataDir)
    
    SumBIN(dataDirL, desDir)
