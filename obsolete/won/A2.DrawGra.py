#!/usr/bin/python

import sys 
from string import *
import os
sys.path.append("/home/won/WON/")
import Util
import Sequtil
import operator
import ChIP


AllGra = Util.GetDirList("Nha/*.gra",DIR=1)
if len(sys.argv)==2: AllGra=Util.GetDirList(sys.argv[1],DIR=1)
for gra in AllGra:
	#if gra.find("romoter")>0: continue
	fp = open(gra,"r")
	length = 0
	count = 0
	for line in fp.xreadlines():
		if line.find(">")>=0: length=count
		count+=1
		if length!=0: break
	fp.close()
	length = length-2
	
	if gra.find("ack")>0: continue
	cmd = 'R -q --no-save --args infile="%s" length="%d"< drawAndFilter.R'%(gra,length)
	print cmd
	os.system(cmd)



