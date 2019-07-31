#!/usr/bin/python

import sys 
from string import *
import os
import os.path
sys.path.append("/home/won/WON/")
import Util
import math
import ChIP
import time,gzip


infile = sys.argv[1]
outdirec = sys.argv[2]
binsize = int(sys.argv[3])


if not os.path.isdir(outdirec):
	os.mkdir(outdirec)

fp = open(infile,"r")
pos = 0
pre_chr = ""
for line in fp.xreadlines():
	try:
		chr, start, stop, val = line.split()
		start = int(start)
		stop = int(stop)
	except: continue
		
	if pre_chr!=chr:
		if pre_chr!="": 
			fout.close()
			sys.exit(1)
		fout = open(outdirec+"/%s.bin"%chr,"w")
		for i in xrange(0,start-binsize+1,binsize):
			fout.write("0\n")
		for j in xrange(start-binsize+1, stop-binsize+1,binsize):
			fout.write(val+"\n")
	for j in xrange(start-binsize+1, stop-binsize+1,binsize):
		fout.write(val+"\n")

fout.close()
	




sys.exit(1)




direc = "hMSC/"
outdirec = "HomerTmp/"

File = Util.GetDirList(direc+"*.bed.gz",DIR=1)

for filename in File:
	fname = filename.split("/")[-1].split(".")[0]
	outfile = filename.replace(".gz","")
	cmd = 'zcat %s | grep -v "track" >%s'%(filename,outfile)
	os.system(cmd)

	tmpdirec = outdirec+fname
	cmdDirec = 'makeTagDirectory %s %s'%(tmpdirec,outfile)
	print cmdDirec
	os.system(cmdDirec)
	cmdBed = 'makeUCSCfile %s  -o %s'%(tmpdirec,fname)
	os.system(cmdBed)

