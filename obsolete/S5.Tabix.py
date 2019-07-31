#!/usr/bin/env python
#$ -S /gpfs/fs0/share/python2.6.5/Python-2.6.5/python

import sys, string, os, glob
from os import path
import pipeUtil
import genomeUtil


def paste_Dir(srcDirL, desDir, chrStrL, nameD, start, step):
    os.system('mkdir -p %s' % desDir)


    fout = open(path.join(desDir, 'header.txt'), 'w')
    for (src, name) in namePairL:
        fout.write('%s\t%s\n' % (src, name)) #nameD[fout.write('\n'.join([path.basename(srcDir) for srcDir in srcDirL])+'\n')
    fout.close()
    
    for chrStr in chrStrL:
        print 'Paste %s' % chrStr
        binFileL = [path.join(binDir, chrStr+'.nrm') for binDir in binDirL ]
        gawkCmd = '\'BEGIN{pos=%d}{OFS="\\t"; print "%s",pos,$0; pos=pos+%d}\'' \
                  % (start, chrStr, step)        
        os.system('paste %s | gawk %s > %s' % \
                  (' '.join(binFileL), gawkCmd, path.join(desDir, chrStr+'.bin')))
                

# TABIX
# cat *.bin | bgzip -c > %s && tabix -s 1 -b 2 -e 2 %s


# Generate Tabix index by pasting given foreground data
#   exclude any data that are used as background
#
# sys.argv[1]: srcDir
# sys.argv[2]: desDir
# sys.argv[3]: tabix name
if __name__ == '__main__':
    if len(sys.argv)<3:
        print 'Arguments: srcDir, desDir, (tabixName)'
        sys.exit(1)                

##    genome, binsize, clonal_mode, fragLen, fgD, namePairL = \
##                pipeUtil.readConfig(LEVEL=1)[0:6]                
    genome, binsize, clonal_mode, fragLen, shiftMode, fgD, namePairL, fragLenD = \
            pipeUtil.readConfig(LEVEL=1)
    
    srcDir = sys.argv[1]
    desDir = sys.argv[2]                
    if len(sys.argv)<4: tbxName = path.basename(os.getcwd())
    else: tbxName = sys.argv[3]
        
    NofChr = genomeUtil.getChrCnt(genome)
    
##    print namePairL
    for (src, name) in sorted(namePairL):
        if src in fgD: namePairL.remove([src, name])
##    binDirL = [src for (src, name) in namepairL]

    # Paste    
    binDirL = [path.join(srcDir, src) for (src, name) in namePairL] #sorted(glob.glob(path.join(srcDir, '*')))
    chrStrL = [genomeUtil.getChrStr(i+1, genome) for i in xrange(NofChr)]
    paste_Dir(binDirL, desDir, chrStrL, namePairL, (binsize+1)/2, binsize)

    # Tabix
    print 'Generating Tabix Index: %s' % tbxName
    tbxPath = path.join(desDir, tbxName+'.gz')
    os.system('cat %s/*.bin | bgzip -c > %s && tabix -s 1 -b 2 -e 2 %s' % \
              (desDir, tbxPath, tbxPath))
    
