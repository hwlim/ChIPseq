#!/usr/bin/env python
#$ -S /gpfs/fs0/share/python2.6.5/Python-2.6.5/python

import sys, string, os
import pipeUtil
import genomeUtil

# sys.argv[1]: genome
# sys.argv[2]: binsize
# sys.argv[3]: srcDir, (./SUB)
# sys.argv[4]: destination directory (./Paste)
# sys.argv[5]: cell
# sys.argv[6:]: source file directories
if __name__ == '__main__':
    genome, binsize = sys.argv[1], int(sys.argv[2])
    NofChr = genomeUtil.getChrCnt(genome)
    srcDir = sys.argv[3]
    pasteDir = sys.argv[4]
    cell = sys.argv[5]
    fgDirL = sys.argv[6:]
    infoDir = './Info'


    # 'cell'.histone: histone list
    os.system('mkdir -p Info')
    hL = sorted([pipeUtil.splitName(fgDir)[1] for fgDir in fgDirL])
    fout = open(os.path.join(infoDir, '%s.histone' % cell), 'w')
    fout.write('\n'.join(hL)+'\n')
    fout.close()

    # cell directory
    os.system('mkdir -p %s' % os.path.join(pasteDir, cell))
    os.system('rm %s' % os.path.join(pasteDir, cell, '*.bin'))
    for i in xrange(1, NofChr): # ignore chrY
        chrStr = genomeUtil.getChrStr(i, genome)
        desFile = os.path.join(pasteDir, cell, chrStr+'.bin')
        fgFileL = [os.path.join(srcDir, fgDir, chrStr+'.nrm') for fgDir in fgDirL]
        cmd = "paste %s | gawk 'BEGIN{pos=1}{OFS=\"\t\"; print pos, $0; pos+=%d;}' > %s" \
              % (' '.join(fgFileL), binsize, desFile)
        print cmd        
        os.system(cmd)
