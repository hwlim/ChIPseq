#!/usr/bin/python

import sys 
from string import *
import os
import os.path
sys.path.append("/home/won/WON/")
import Util
import math
import ChIP
import time
import operator

binsize=20
#direc = "Homer/"
direc = "HomerTmp/"
outdirec = "BIN/"

All = Util.GetDirList(direc+"*",DIR=1)
#All = [x for x in All if x.find("Breast_Myoepithelial_")>0]
Cell = [x.split("/")[-1] for x in All if x.find("ccdc")>=0]
Cell = list(set(Cell))
Cell.sort()
print Cell

fout = open("name.cell","w")
fout.write("\n".join(Cell))
fout.close()


T = []
for ch in xrange(1,25):
	chr = "chr"+str(ch)
	if ch==23: chr="chrX"
	elif ch==24: chr="chrY"
	
	MergeList = []
	outfile = outdirec+"%s.bin"%chr
	for cell in Cell:
		fhisfile = direc+"%s/%s.tags.nrm"%(cell,chr)
		MergeList.append(fhisfile)
	print MergeList

	cmdmerge = "paste "
	for merge in MergeList:
		cmdmerge += merge+" "
	cmdmerge += "| gawk '{s+=%d} {print "%(binsize)
	cmdmerge += '"%s\\t"s-%d"\\t"'%(chr,binsize/2)
	for m in range(1,len(MergeList)+1):
		cmdmerge += '$%d"\\t"'%m
	cmdmerge += "}' > %s"%outfile
	print cmdmerge
	os.system(cmdmerge)
	sys.exit(1)

