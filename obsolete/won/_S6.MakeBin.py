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

F = Util.GetDirList("HomerTmp/*",DIR=1)
F = [x for x in F if (os.path.isdir(x)) & (x.find("119_2")<0)] # remove background
shortF = [x.split("/")[-1] for x in F]

fout = open("name.cell","w")
fout.write("\n".join(shortF))
fout.close()


fout = open("HomerTmp/hMSC.bin","w")
for i in xrange(1,25):
	if i==23: chr="chrX"
	elif i==24: chr="chrY"
	else: chr="chr"+str(i)
	print chr	
	fname = "%s.tags.nrm"%(chr)
	T = []
	for j in xrange(len(F)):
		name = F[j]+"/"+fname
		L = Util.ReadColumn(name,[0])
		T.append(L)
		L = []
	O = zip(*T)
	T = []


	for c in xrange(len(O)):
		pnt = c*20+10
		dO = [chr, pnt] + list(O[c])
		dO = map(str,dO)

		fout.write("\t".join(dO)+"\n")
	O = []
fout.close()




