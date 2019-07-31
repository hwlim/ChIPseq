#!/usr/bin/env python

import sys, string, os, glob
from os import path
import pipeUtil
from optparse import OptionParser
from multiprocessing import Pool
##from operator import itemgetter


def HomerPeak(optD):
    try:
        optStr=''
##        print optD['name']
        if optD.has_key('name'): optStr = '-n %s' % optD['name']
        if optD.has_key('outDir'):
            optStr = optStr + ' -o %s' % path.normpath(optD['outDir'])
        else:
            optStr = optStr + ' -o HomerPeak/%s' \
                     % path.basename(path.normpath(optD['chip']))
            
        if optD.has_key('optionStr'): optStr = optStr + ' -e "%s"' % optD['optionStr']
        
        cmd = 'runHomerPeak.sh %s %s %s' % \
              (optStr, path.normpath(optD['chip']), path.normpath(optD['input']))
        print cmd
        os.system(cmd)
        return 'Success'
    except KeyboardInterrupt, e:
        pass





# Command line argument: Number of process
if __name__ == '__main__':
    parser=OptionParser(description='Wrapper program for Homer findPeaks',\
                        usage='%prog [options]')
    parser.add_option('-c',  dest='confFile', default='HomerPeakConfig.txt',
                      help='Configuration file')
    parser.add_option('-p',  dest='MaxProc', default=1,
                      help='Maximum number of processes.')

    if len(sys.argv)==1:
        parser.print_help()
        exit(-1)
        
    (opts, args) = parser.parse_args()

    for bed in args:
        if not os.path.isfile(confFile):
            print >> sys.stderr, "Error: %s does not exist." % bed
            exit(-1)
    

    MaxProc = int(opts.MaxProc)
    configL = pipeUtil.readHomerPeakConfig(opts.confFile)

    for config in configL:
        if not os.path.isdir(config['chip']):
            print >> sys.stderr, "Error: Tag Dir %s does not exist." % config['chip']
            exit(-1)
        if not config['input']!='none' and os.path.isdir(config['input']):
            print >> sys.stderr, "Error: Input Dir %s does not exist." % config['input']
            exit(-1)

    print 'Starting multiprocessing: %d' % MaxProc
    pool = Pool(processes=MaxProc)
    p = pool.map_async(HomerPeak, configL)
    try:
        results = p.get(0xFFFF)
    except KeyboardInterrupt:
        print 'Keyboard Interrupt by Ctrl-c, terminating workers...'
        pool.terminate()
        pool.join()
    else:
        print 'Quitting normally'
        pool.close()
        pool.join()

