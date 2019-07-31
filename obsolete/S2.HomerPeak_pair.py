#!/usr/bin/env python

import sys, string, os, glob
from os import path
import pipeUtil
from optparse import OptionParser
from multiprocessing import Pool
##from operator import itemgetter


def HomerPeak((chipDir, inputDir, name, fdr)):
    try:
        cmd = 'run_homer.sh %s %s %s %f' % (chipDir, inputDir, name, fdr)
##        print cmd
        os.system(cmd)
        return 'Success'
    except KeyboardInterrupt, e:
        pass





# Command line argument: Number of process
if __name__ == '__main__':
    if len(sys.argv)<4:
        print 'Usage: S2.HomerPeak_pair.py [MaxProc] [desDir] [tagDir]'
        print '\tPairwise differential peak finding'
        sys.exit(1)

    MaxProc = int(sys.argv[1])
    desDir = sys.argv[2]
    tagDir = path.normpath(path.join(os.getcwd(), sys.argv[3]))

    commonOptD, dataOptD = pipeUtil.readConfig2()

    pwd=os.getcwd()
    os.system('mkdir -p %s' % desDir)
    print 'Changing directory: %s' % desDir
    os.chdir('%s' % desDir)

##    print fgD

    fgL = [ fg for fg in sorted(dataOptD) ]
    argL = []
    for i, fg  in enumerate(fgL[:-1]):
        for j, bg in enumerate(fgL[i+1:]):
            argL = argL + \
                   [(path.join('..', tagDir, fg), \
                     path.join('..', tagDir, bg), \
                     dataOptD[fg]['name']+'_Over_'+dataOptD[bg]['name'],\
                     0.001), \
                    (path.join('..', tagDir, bg), \
                     path.join('..', tagDir, fg), \
                     dataOptD[bg]['name']+'_Over_'+dataOptD[fg]['name'],\
                     0.001) \
                    ]
    
##    argL = [(path.join('..', tagDir, fg), \
##             path.join('..', tagDir, dataOptD[fg]['bg']) if dataOptD[fg]['bg'].lower() != "none" else "none",\
##             dataOptD[fg]['name'],\
##             0.001) \
##            for fg in sorted(dataOptD)]

##    for arg in argL:
##        print arg
    print 'Starting multiprocessing: %d' % MaxProc
    pool = Pool(processes=MaxProc)
    p = pool.map_async(HomerPeak, argL)
    try:
        results = p.get(0xFFFF)
    except KeyboardInterrupt:
        print 'Keyboard Interrupt by Ctrl-c'
        os.chdir('..')

##        pool.terminate()
    print 'Changing directory: %s' % pwd
    os.chdir('..')
