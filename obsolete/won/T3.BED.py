#!/usr/bin/python

import sys 
from string import *
import os
sys.path.append("/home/won/WON/")
import Util
import Sequtil
import operator

wide = 100
maxnum = 1000
outdirec = "Nha/"
F = Util.GetDirList("/home/imsi/ccdc/*.txt",DIR=1)
for fname in F:
	if fname.find("BiModel")<0:
		Input = Util.ReadColumn(fname,[0,1,2,3])
	else:
		Input = Util.ReadColumn(fname,[0,4,3,2])
	if Input[0][0].find("chr")<0: Input=Input[1:]
	Input.sort(lambda x,y:cmp(float(y[2]),float(x[2])))

	count = 0
	
	filename = fname.split("/")[-1]
	fout = open(outdirec+"%s.bed"%filename,"w")
	count = 0
	for input in Input:
		chr, pos,intensity, intensity2 = input
		if chr=="chr23": chr="chrX"
		if chr=="chr24": chr="chrY"
		pos = int(pos)
		O = [chr, pos-wide, pos+wide]
		O = map(str,O)
		fout.write("\t".join(O)+"\n")
		count += 1
		if count == maxnum: break

	fout.close()
