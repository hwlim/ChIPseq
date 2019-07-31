#!/usr/bin/env python
import os, sys, glob, string
##from multiprocessing import Pool
import pipeUtil
import multiprocessing

# NEED TO CHANGE TO MULTIPROCESSING VERSION
if __name__=='__main__':

    MQScrTh = 30
    srcPath = './BAM/'  
    desPath = './BED/'

    bamfile = sys.argv[1].split('/')[-1]
  
    bamDataL = glob.glob(srcPath+'*.bam')
    bamDataL = [bamData.replace('\\', '/') for bamData in bamDataL] # window compatibility

    cmd = 'mkdir -p BED'
    os.system(cmd)

    try:
        pool.map(BAM2BED, bamDataL)
    except KeyboardInterrupt:
        pool.terminate()
