#!/usr/bin/env python
import sys, string, os
##from string import *
from os import path
from multiprocessing import Pool
import pipeUtil
import FGvsBG

##sys.path.append("/home/heewlim/Research/PythonModules/Won/")
##import Util
##import math
##import ChIP
##import time
##import operator




####def GetHistoneList(direc,cell):
####	All = Util.GetDirList(direc+"*.bin",DIR=1)
####	All = [x for x in All if cell in x.split(".")]
####
####	Histone = [x.split(".")[2] for x in All]
####	Histone = [x for x in Histone if x.find("H")==0]
####	Histone = list(set(Histone))
####	Histone = ChIP.HistoneSort(Histone)  # arrange histone name
####
####	Dnase = [x for x in All if "ChromatinAccessibility" in x.split(".")]
####	if len(Dnase)>0: Histone.append("ChromatinAccessibility")
####	print "All Histone in %s:"%cell,Histone
####	return Histone
##
##def Sub(H,I):
##	O = H[:]
##	for i in xrange(len(O)):
##		O[i] = float(H[i])-float(I[i])
##		#if O[i]<0: O[i]=0  # from Nha's suggestion 
##	return O
##
##def Div(O,I):
##	O = H[:]
##	for i in xrange(len(H)):
##		O[i] = (float(H[i])+5.)/(float(I[i])+5.) -1.0
##		if O[i]<0: O[i]=0
##	return O
##
##def Average(Total):
##	denom = float(len(Total))
##	if denom==1: return Total[0]  
##	O = []
##	for i in xrange(len(Total[0])):
##		avg = sum(map(operator.itemgetter(i),Total))/denom
##		O.append(avg)
##	return O
##
##def GetInfo(fname):
##	fp = open(fname,"r")
##	line = fp.readline()
##	if line.split()[0]!="name": 
##		fp.close(); print "error 0"
##		return 0,0
##	try:
##		line1 = fp.readline()
##		totaltag = float(line1.split()[2])
##		line2 = fp.readline()
##		fragment = int(line2.split("=")[1])
##	except: 
##		fp.close()
##		print "error 1"
##		return 0,0
##
##	return totaltag, fragment
##
##
##def RPKM(T,total):
##	factor = total/1000000. #factor = 2M counts
##	O = T[:]
##	for j in xrange(len(T)):
##		O[j] = float(T[j])/factor/(binsize/1000.)
##	return O
##
##def PrintHistone(T,fname):
##	print "writing", fname
##	fp = open(fname,"w")
##	for i in xrange(len(T)):
##		fp.write("%f\n"%(T[i]))
##	fp.close()
##	return
	

def Job((fgL, bg, srcDir, desDir)):
    print 'FGvsBG: Subtracting background %s' % bg
    FGvsBG.FGvsBG(fgL, bg, srcDir, desDir)

def Job2((fgOnly, srcDir, desDir)):
    print 'FGvgBG: Copying %s/%s --> %s/' % (srcDir, fgOnly, desDir)
    os.system('cp -r %s/%s %s/%s' % (srcDir, fgOnly, desDir, fgOnly))
    
def main():
    if len(sys.argv)<4:
        print('Usage:\n\tS4.FGvsBG [MaxProc] [desDir] [srcDir]')
        sys.exit(1)

    maxProc = int(sys.argv[1])
    desDir = sys.argv[2]
    srcDir = sys.argv[3]

    ##genome, binsize, MaxProc, fgL, bgL = pipeUtil.readConfig()
    genome, binsize, clonal_mode, fragLen, shiftMode, fgD, namePairL, fragLenD = \
            pipeUtil.readConfig(LEVEL=1)

    if len(fgD)==0:
        print 'No data.\nNothing to be done.'
        sys.exit(1)

    os.system('mkdir -p %s' % desDir)


    pool = Pool(processes=maxProc)

    # copy data sets that has no background data
    fgOnly = fgD[None]
    for fg in sorted(fgOnly):
        if fg in fgD.keys(): fgOnly.remove(fg)
    
    argL = [[fg, srcDir, desDir] for fg in sorted(fgD[None])]
    try: pool.map(Job2, argL)
    except KeyboardInterrupt: pool.terminiate()
    
    # back ground substraction
    del fgbgD[None]
    argL = [[fgD[bg], bg, srcDir, desDir] for bg in sorted(fgD)]
    try: pool.map(Job, argL)
    except KeyboardInterrupt: pool.terminiate()




##fgD = {}
##for i, bg in enumerate(bgL):
##    if not fgD.has_key(bg): fgD[bg] = []
##    fgD[bg].append(fgL[i])
##
####        bgD = dict([fgL[i], bgL[i]] for i in xrange(len(fgL)))
##
##for bg in bgD:
##    bgAll = Util.GetDirList(nrmDir+bg+"/*.nrm",DIR=1)
##    fgL = fgD[bg]
##for bgChr in bgAll:
##    nrmBG = Util.ReadColumn(bgChr,[0])
##for fg in fgL:
##    fgChr = bgChr.replace(bg, fg)
##    nrmFG = Util.ReadColumn(fgChr, [0])
##
##subName = fgChr.replace(srcDir, desDir)
####divName = subName.replace(".nrm", ".div")
##
##LevelH = Sub(nrmFG,nrmBG)
####DivH = Div(nrmFG,nrmBG)
##PrintHistone(LevelH,subName)
##PrintHistone(DivH,divName)
                                


if __name__ == '__main__':
        main()
