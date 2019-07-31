#!/usr/bin/env python
#$ -S /gpfs/fs0/share/python2.6.5/Python-2.6.5/python

import sys, string, os, glob
from os import path
##import pipeUtil
import genomeUtil


def pasteBIN_Dir(srcDirL, desDir, infoFile, chrStrL, start, step):
    os.system('mkdir -p %s' % desDir)
    
    fout = open(infoFile, 'w')
    fout.write('\n'.join([path.basename(srcDir) for srcDir in srcDirL])+'\n')
    fout.close()

    for chrStr in chrStrL:
        print 'Paste %s' % chrStr
        binFileL = [path.join(binDir, chrStr+'.for.bin') for binDir in binDirL ]
        gawkCmd = '\'BEGIN{pos=%d}{OFS="\\t"; print "%s",pos,$0; pos=pos+%d}\'' \
                  % (start, chrStr, step)        
        os.system('paste %s | gawk %s > %s' % \
                  (' '.join(binFileL), gawkCmd, path.join(desDir, chrStr+'.for.bin')))
                
        binFileL = [path.join(binDir, chrStr+'.rev.bin') for binDir in binDirL ]
        gawkCmd = '\'BEGIN{pos=%d}{OFS="\\t"; print "%s",pos,$0; pos=pos+%d}\'' \
                  % (start, chrStr, step)
        os.system('paste %s | gawk %s > %s' % \
                  (' '.join(binFileL), gawkCmd, path.join(desDir, chrStr+'.rev.bin')))
    

# TABIX
# cat *.bin | bgzip -c %s && tabix -s 1 -b 2 -e 2 %s

# sys.argv[1]: genome
# sys.argv[2]: binsize
# sys.argv[3]: srcDir
# sys.argv[4]: desDir
if __name__ == '__main__':
    if len(sys.argv)<5:
        print 'Arguments: genome, binsize, srcDir, desDir'
        sys.exit(1)

    genome, binsize = sys.argv[1], int(sys.argv[2])
    srcDir = sys.argv[3]
    desDir = sys.argv[4]

    NofChr = genomeUtil.getChrCnt(genome)

    
    binDirL = sorted(glob.glob(path.join(srcDir, '*')))
    chrStrL = [genomeUtil.getChrStr(i+1, genome) for i in xrange(NofChr)]
    pasteBIN_Dir(binDirL, desDir, path.join(desDir, 'info.txt'), \
                 chrStrL, binsize/2, binsize)
