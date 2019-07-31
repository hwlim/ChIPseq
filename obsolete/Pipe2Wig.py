#!/usr/bin/env python
# Extract intermediate results to a destination directory for debugging purpose

import glob, os, sys, string
from os import path
import util
import genomeUtil
import pipeUtil

binDir = './BIN'
nrmDir = './NRM'
subDir = './SUB'
desDir = './test'

def Wigfy(srcFile, desDir, dataName, chrStr, binsize, option):

    if option=='BIN': postfix = 'bin'
    elif option=='NRM': postfix = 'nrm'
    elif option=='SUB': postfix = 'sub'
    else: raise RuntimeError('Invalid option: %s' % option)
    
    desFile = path.join(desDir, '.'.join([dataName, postfix, chrStr, 'wig']))
    wigName = dataName + '.' + postfix
    print srcFile, '-->', desFile
    genomeUtil.writeWigHeader(desFile, chrStr, 1, binsize, None, wigName, wigName)
    os.system('cat %s >> %s' % (srcFile, desFile))


genome, binsize, clonal_mode, fgD = pipeUtil.readConfigSGE(LEVEL=1)
fg = sys.argv[1]

bgName = None
for bg in fgD:
    if fg in fgD[bg]: 
        bgName = bg
        break
if not bgName:
    raise RuntimeError('No such FG data in config.txt: %s' % fg)

chrStr = 'chr1'
fgSrcFile = path.join(binDir, fg, 'chr1.bin')
Wigfy(fgSrcFile, desDir, fg, chrStr, binsize, 'BIN')
bgSrcFile = path.join(binDir, bg, 'chr1.bin')
Wigfy(bgSrcFile, desDir, bg, chrStr, binsize, 'BIN')

fgSrcFile = path.join(nrmDir, fg, 'chr1.nrm')
Wigfy(fgSrcFile, desDir, fg, chrStr, binsize, 'NRM')
bgSrcFile = path.join(nrmDir, bg, 'chr1.nrm')
Wigfy(bgSrcFile, desDir, bg, chrStr, binsize, 'NRM')

c, h, repN = pipeUtil.splitName(fg)
srcFile = path.join(subDir, c, h, 'chr1.nrm')
Wigfy(srcFile, desDir, fg, chrStr, binsize, 'SUB')
