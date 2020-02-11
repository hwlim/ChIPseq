#!/usr/bin/python

import sys 
from string import *
import os
sys.path.append("/home/won/WON/")
import Util
import Sequtil
import operator
import multiprocessing

Fa = Util.GetDirList("15*.fa",DIR=1)
count =0
MaxProc=6


def Job(fa):
	out = "%s.meme"%(fa)
	cmd = "~/bin/meme %s -dna -minw 6 -maxw 15 -text -revcomp -nmotifs 9 -mod oops -maxsize 1000000 > %s"%(fa,out)
	print cmd
	os.system(cmd)
	return

for fa in Fa:
	#if fa.find("type")<0: continue
	
	p = multiprocessing.Process(target=Job,args=(fa,))
	p.start()
	if count%MaxProc==(MaxProc-1): p.join()
	print count,"/", len(Fa)
	count +=1

