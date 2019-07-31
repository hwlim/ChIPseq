#!/usr/bin/env python
#$ -S /gpfs/fs0/share/python2.6.5/Python-2.6.5/python
# GRO-seq BED tag coverage check
# ASSUMPTION: bed file is already sorted by chr, direction (optionaly start position)
# Special check script for MITCH data
import sys, os, string, glob
from os import path

def gsCheckBEDCvg(srcFile):
    outFile = path.join(path.dirname(srcFile),\
                        path.basename(srcFile).replace('.bed','')+'.info')
    print 'Check BED Coverage: %s' % srcFile
    bedCoverD={}
    fin = open(srcFile, 'r')
    for bed in fin:
        field = bed.split()
        if len(field)==0 or field[0][0:3]!='chr': continue
        chrStr, start, end, direc = \
                field[0], int(field[1]), int(field[2]), field[5]
        length = end-start+1
        if not bedCoverD.has_key(length): bedCoverD[length] = {}
        if not bedCoverD[length].has_key(chrStr): bedCoverD[length][chrStr] = []
        if direc=='+': bedCoverD[length][chrStr].append(start)
        elif direc=='-': bedCoverD[length][chrStr].append(-end)
        else: raise RuntimeError('Invalid direction: %s' % direc)
    fin.close()
    
    coverPosSD = {}
    preN = 0
    fout = open(outFile, 'w')
    fout.write('ReadLength\tTagCnt\tAccuPosCnt\tIncrease\tIncreaseRatio\n')
    print 'ReadLength\tTagCnt\tAccuPosCnt\tIncrease\tIncreaseRatio'
    for i, length in enumerate(sorted(bedCoverD, reverse=True)):
        curN = 0
        tagN = 0
        for chrStr in bedCoverD[length]:
            if not coverPosSD.has_key(chrStr): coverPosSD[chrStr] = set()
            coverPosSD[chrStr] = coverPosSD[chrStr].union(set(bedCoverD[length][chrStr]))
            tagN = tagN + len(set(bedCoverD[length][chrStr]))
            curN = curN + len(coverPosSD[chrStr])
        fout.write('%d\t%d\t%d\t%d\t%.3f\n' %(length, tagN, curN, curN-preN, float(curN-preN)/curN))
        print '%d\t%d\t%d\t%d\t%.3f' %(length, tagN, curN, curN-preN, float(curN-preN)/curN)
        preN = curN
    fout.close()


def gsCheckBEDCvg_Dir(srcDir):
    fileL = glob.glob(path.join(srcDir, '*.bed'))
    for bedFile in sorted(fileL):
##        dataName = path.basename(bedFile).replace('.bed', '')
        gsCheckBEDCvg(bedFile)
    
        
if __name__=='__main__':
    
    if len(sys.argv)<2:
        print 'Arguments: src'
        sys.exit(1)
    src = sys.argv[1]

    if path.isdir(src):
        gsCheckBEDCvg_Dir(src)
    else:
        gsCheckBEDCvg(src)
