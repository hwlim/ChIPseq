#!/usr/bin/env python

import sys, string, os, glob
from optparse import OptionParser
from os import path
import pipeUtil
from multiprocessing import Pool


def BED2TSV((srcFile, desDir, genome)):
    try:
        dataName = path.basename(srcFile).replace('.bed', '')
        tsvDir = path.join(desDir, dataName)
        cmd = 'makeTagDirectory %s %s -format bed -force5th -genome %s -checkGC 2> %s' % \
              (tsvDir, srcFile, genome, tsvDir+'.log')
        print >> sys.stderr, cmd
        os.system(cmd)
        cmd = 'drawGCplot.sh %s' % (tsvDir)
        print >> sys.stderr, cmd
        os.system(cmd)
        cmd = 'drawTCDistrib.sh %s' % (tsvDir)
        print >> sys.stderr, cmd
        os.system(cmd)
        cmd = 'drawTLDistrib.sh %s' % (tsvDir)
        print >> sys.stderr, cmd
        os.system(cmd)
        return 'Success'
    except KeyboardInterrupt, e:
        pass



# Command line argument: Number of process
if __name__ == '__main__':
    parser=OptionParser(description='Wrapper program for Homer\'s makeTagDirectory for batch processing.',\
                        usage='%prog [options] [bedFile1] (bedFile2) ...')
    parser.add_option('-g',  dest='genome', choices=['mm8','mm9','hg18','hg19'],
                      help='Target genome, mm8/mm9/hg18/hg19')
    parser.add_option('-d',  dest='desDir', default='TSV',
                      help='Destination directory under which a tag directory is created (for batch mode).')
##    parser.add_option('-o',  dest='tsvDir', default='TSV',
##                      help='Tag directory (for single mode).')
    parser.add_option('-p',  dest='MaxProc', default=1,
                      help='Maximum number of processes.')

    (opts, args) = parser.parse_args()

    for m in ['genome']:
        if not opts.__dict__[m]:
            print >> sys.stderr, "Mandatory option is missing: %s" % m
            parser.print_help()
            exit(-1)
    if len(args)==0:
        print >> sys.stderr, "bed2tsv.py requires at least one input bed file."
        exit(-1)

    for bed in args:
        if not os.path.isfile(bed):
            print >> sys.stderr, "Error: %s does not exist." % bed
            exit(-1)
            
    os.system('mkdir -p %s' % opts.desDir)

    argL = [[bedFile, opts.desDir, opts.genome] for bedFile in args]

    print >> sys.stderr, 'Starting multiprocessing: %d' % int(opts.MaxProc)
    pool = Pool(processes=int(opts.MaxProc))
    p = pool.map_async(BED2TSV, argL)
    try:
        results = p.get(0xFFFF)
    except KeyboardInterrupt:
        print >> sys.stderr, 'Keyboard Interrupt by Ctrl-c'
