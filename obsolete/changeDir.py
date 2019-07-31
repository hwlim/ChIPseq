#!/usr/bin/env python
# reorganize directory name and structure
# flat name --> cell/histone hierarchy, normalized histone name
import glob, os, sys, string
import pipeUtil

srcDir = './SUB'
dirL = glob.glob('./SUB/*')

for d in dirL:
    dataName = os.path.basename(d)
    c, h, repN = pipeUtil.splitName(dataName)
    desDir = os.path.join(srcDir, c, h)
    print desDir
    os.system('mkdir -p %s' % os.path.join(srcDir, c, h))
    os.system('mv %s %s'% (os.path.join(d, '*.nrm'), desDir))
    os.system('rmdir %s'% d)
