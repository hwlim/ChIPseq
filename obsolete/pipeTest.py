#!/usr/bin/env python

##import sys
##sys.path.append('/home/heewlim/Research/PythonModules/')


##from multiprocessing import Pool
##import os
##
##def f(x):
##    print '< %d, %s, %d >' % (os.getpid(), x[0], x[1])
##    
##
##if __name__=='__main__':
##    pool = Pool(processes=3)
##
##    myArg = [['a',1], ['b',2], ['c',3], ['d',4], ['e',5], ['f',6], ['g',7]]
##    pool.map(f, myArg);

import sys
import genomeUtil
import pipeUtil

print pipeUtil.readConfigSGE()
##
##import pipeUtil
##
##print pipeUtil.readHomerInfo('/ifs/h/heewlim/Pipeline/ENCODE/HomerTmp/wgEncodeBroadHistoneK562H3k9me3StdAlnRep2/')
