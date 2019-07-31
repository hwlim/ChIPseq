#!/usr/bin/env python
#$ -S /gpfs/fs0/share/python2.6.5/Python-2.6.5/python
import sys, os, string, glob
import pipeUtil
import util



# RPKM
def BIN2NRM(srcFile, desFile, binsize, totalTagCnt):
    factor1 = totalTagCnt/1000000. #factor = 2M counts
    factor2 = binsize/1000.
    factor = factor1*factor2
    
##    print "writing %s" % desFile
    fin = open(srcFile, 'r')
    fout = open(desFile, 'w')
    for line in fin:
        fout.write('%.3f\n' % (float(line)/factor))

    fin.close()
    fout.close()

def BIN2NRM_DIR(srcDir, desDir, binsize):
    print 'BIN2NRM: %s --> %s' % (srcDir, desDir)
    os.system('mkdir -p %s' % desDir)
    totalTagCnt = pipeUtil.readBINInfo(srcDir)
    binfileL = glob.glob(os.path.join(srcDir,'*.bin'))
    for binfile in sorted(binfileL):
        nrmName=os.path.basename(binfile).replace('.bin','.nrm')
        BIN2NRM(binfile, os.path.join(desDir,nrmName), binsize, totalTagCnt)


# command line arguments:
# sys.argv[1]: srcDir or srcFile
# sys.argv[2]: desDir or srcFile
if __name__=='__main__':
    if len(sys.argv)<3:
        print 'Usage:\n\tBIN2NRM.py srcDir desDir'
        sys.exit(1)
    
    genome, binsize, clonal_mode, fragLen = pipeUtil.readConfig(LEVEL=0)
    
    src, des = sys.argv[1], sys.argv[2]
    if os.path.isdir(src):
        BIN2NRM_DIR(src, des, binsize)
    else:
        BIN2NRM(src, des, binsize)
    
