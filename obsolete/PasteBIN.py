#!/usr/bin/env python
#$ -S /gpfs/fs0/share/python2.6.5/Python-2.6.5/python

import sys, string, os, glob
from os import path
##import pipeUtil
import genomeUtil


def pasteBIN_Dir(srcDirL, desDir, chrStrL, nameD, start, step):
    os.system('mkdir -p %s' % desDir)


    fout = open(path.join(desDir, 'header.txt', 'w'))
    for src in nameD:
        fout.write('%s\t%s\n' % (src, nameD[src])) #nameD[fout.write('\n'.join([path.basename(srcDir) for srcDir in srcDirL])+'\n')
    fout.close()
    
##    for chrStr in chrStrL:
##        print 'Paste %s' % chrStr
##        binFileL = [path.join(binDir, chrStr+'.bin') for binDir in binDirL ]
##        gawkCmd = '\'BEGIN{pos=%d}{OFS="\\t"; print "%s",pos,$0; pos=pos+%d}\'' \
##                  % (start, chrStr, step)        
##        os.system('paste %s | gawk %s > %s' % \
##                  (' '.join(binFileL), gawkCmd, path.join(desDir, chrStr+'.bin')))
                
   

# TABIX
# cat *.bin | bgzip -c %s && tabix -s 1 -b 2 -e 2 %s

#### sys.argv[1]: genome
#### sys.argv[2]: binsize
# sys.argv[1]: srcDir
# sys.argv[2]: desDir
if __name__ == '__main__':
    if len(sys.argv)<3:
        print 'Arguments: srcDir, desDir'
        sys.exit(1)                

    genome, binsize, clonal_mode, fragLen, fgD, nameD = \
                pipeUtil.readConfig(LEVEL=1)[0:6]                
    srcDir = sys.argv[3]
    desDir = sys.argv[4]                

    NofChr = genomeUtil.getChrCnt(genome)

    binDirL = nameD.keys()
    for src in sorted(nameD):
        if src in fgD: del nameD[src]
    
    binDirL = [path.join(srcDir, src) for src in nameD] #sorted(glob.glob(path.join(srcDir, '*')))
    chrStrL = [genomeUtil.getChrStr(i+1, genome) for i in xrange(NofChr)]
    pasteBIN_Dir(binDirL, desDir, chrStrL, nameD, (binsize+1)/2, binsize)
