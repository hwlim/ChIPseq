#!/usr/bin/env python
#$ -S /gpfs/fs0/share/python2.6.5/Python-2.6.5/python
#
# TSV --> BIN file

import sys, os, string, glob, math
from os import path
import pipeUtil
import genomeUtil
import Util

# overlap length check, including boundary
def checkOver(bed, tag):
    if bed[1] <= tag[0]:   # bed --> tag, no overlap
        flag = -2
        overLen = 0   
    elif bed[0] >= tag[1]: # tag --> bed, no overlap
        flag = -1
        overLen = 0        
    else:
        flag = 1
        overLen = min(bed[1]-tag[0],\
                      tag[1]-bed[0],\
                      bed[1]-bed[0],\
                      tag[1]-tag[0])    # including boundary
    return flag, overLen


# how to compute bin average?
# current: overlap length / bin size
# but what if fragment length < bin size?
def BED2BIN(tagL, chrLen, binsize, desFile, avgMode=True):
    fout = open(desFile, 'w')
    curTagIdx = 0
    tagIdxS = set()  # set of tag index, overlapped once
    for binStart in xrange(0, chrLen-binsize, binsize): # binEnd = binStart + binsize -1
        tCnt = 0.0
        for tagIdx in list(tagIdxS):
            flag, overLen = checkOver([binStart, binStart+binsize], tagL[tagIdx][0:2])
            if flag == 1:
                if avgMode: tCnt = tCnt + overLen*tagL[tagIdx][2]
                else: tCnt = tCnt + tagL[tagIdx][2]
            elif flag == -1: tagIdxS.remove(tagIdx)
            elif flag == -2: continue
            else: raise RuntimeError('Invalid flag: %d' % flag)

        while curTagIdx < len(tagL):
            flag, overLen = checkOver([binStart, binStart+binsize], tagL[curTagIdx])
            if flag == -2: break
            elif flag == -1:
                curTagIdx = curTagIdx + 1;
                break
            elif flag == 1:
                if avgMode: tCnt = tCnt + overLen*tagL[curTagIdx][2]
                else: tCnt = tCnt + tagL[curTagIdx][2]
                tagIdxS.add(curTagIdx)
                curTagIdx = curTagIdx + 1
            else: raise RuntimeError('Invalid flag: %d' % flag)
            
        if avgMode: tCnt = tCnt/binsize
        fout.write('%.3f\n' % tCnt)
    fout.close()
    

    
def TSV2BIN(srcFile, desFile, maxCloneCnt, binsize, chrLen, fragLen, shiftMode):
   
    Item = Util.ReadColumn(srcFile,[1,2,3]) #chr, pos, strand, Numtag(*)
    tagCnt = 0.0
    N = []

##    if not shiftMode: shiftMode = 0;
##    print 'shiftMode = %d' % shiftMode
    
    for i in xrange(len(Item)):
        cnt = min(maxCloneCnt, float(Item[i][2])) # determined by clonal_mode

        if shiftMode:
            if Item[i][1]=="0":
                O=[ int(Item[i][0])+fragLen/2-75, int(Item[i][0])+fragLen/2+75, cnt ]
            elif Item[i][1]=="1":
                O=[ int(Item[i][0])-fragLen/2-75, int(Item[i][0])-fragLen/2+75, cnt ]
            else: print "so what"; sys.exit(1)
        else:
            if Item[i][1]=="0": O=[ int(Item[i][0]), int(Item[i][0])+fragLen, cnt ]
            elif Item[i][1]=="1": O=[ int(Item[i][0])-fragLen, int(Item[i][0]), cnt ]
            else: print "so what"; sys.exit(1)
       
        tagCnt = tagCnt + float(cnt)
        N.append(O)    

##    if clonal_mode == 0:    # allow monoclonal reads
##        for i in xrange(len(Item)):
##            if Item[i][1]=="0": O=[ int(Item[i][0]), int(Item[i][0])+fragLen, float(Item[i][2]) ]
##            elif Item[i][1]=="1": O=[ int(Item[i][0])-fragLen, int(Item[i][0]),float(Item[i][2]) ]
##            else: print "so what"; sys.exit(1)
##            tagCnt = tagCnt + int(float(Item[i][2]))
##            N.append(O)
##    elif clonal_mode == 1:  # not allow monoclonal reads
##        for i in xrange(len(Item)):
##            if Item[i][1]=="0": O=[ int(Item[i][0]), int(Item[i][0])+fragLen, 1.0]
##            elif Item[i][1]=="1": O=[ int(Item[i][0])-fragLen, int(Item[i][0]), 1.0]
##            else: print "so what"; sys.exit(1)
##            tagCnt = tagCnt + 1
##            N.append(O)
##
##    totalTagCnt = TSV2BIN.TSV2BIN_DIR(srcDir, desDir, clonal_mode, binsize, chrLenD, fragLen, shiftMode)

    # bed --> bin, file writing
    N.sort(key=lambda x:x[0])
    BED2BIN(N, chrLen, binsize, desFile, avgMode=True)
    
    return tagCnt


# read clonal count information from TSV/tagInfo.txt
def readClonalCntL(srcDir):
    fin = open(path.join(srcDir, 'tagInfo.txt'), 'r')
    for i in xrange(8): fin.readline() # skip unnecessary lines
    cntL = []
    readCntL = []
    for line in fin:
        field = line.split()
        if len(field) != 2: continue
        cntL.append(float(field[0]))
        readCntL.append(float(field[1]))
    fin.close()
    return cntL, readCntL

# calculate maximum allowed clonal count
# current: mean + 3*std
def getMaxClonalCnt(cntL, readCntL):
    m = (sum([cnt*readCnt for cnt, readCnt in zip(cntL, readCntL)])\
        /sum([readCnt for readCnt in readCntL]))
##    print m
    sd = math.sqrt(sum([(cnt-m)**2*readCnt for cnt, readCnt in zip(cntL, readCntL)])\
        /sum([readCnt for readCnt in readCntL]))
    return int(m+3*sd)


	
def TSV2BIN_DIR(srcDir, desDir, clonal_mode, binsize, chrLenD, fragLen=None, shiftMode=None):
    print 'TSV2BIN: %s --> %s' % (srcDir, desDir)
    os.system('mkdir -p %s' % desDir)
    
    if not fragLen:
        total, fragLen, avgTPP = pipeUtil.readTSVInfo(srcDir)
    maxClonalCnt = None
    if clonal_mode==0: maxClonalCnt = float('inf')   # allow monoclonal
    elif clonal_mode >= 1: maxClonalCnt = clonal_mode # not allow monoclonal
##    elif clonal_mode==2:                            # partial allow monoclonal
##        cntL, readCntL = readClonalCntL(srcDir)
##        maxClonalCnt = getMaxClonalCnt(cntL, readCntL)
##        print('clonal_mode=2\nMaxClonalCount=%d' % maxClonalCnt)
##    elif clonal_mode  2:                           # fixed max clonal cnt
##        maxClonalCnt = clonal_mode
    else: raise RuntimeError('Invalid clonal_mode: %d' % clonal_mode)
    
    totalTagCnt = 0.0
    tsvL = glob.glob(path.join(srcDir,'*.tsv'))
    print 'TSV2BIN:fragLen=%d'%fragLen
    for tsv in sorted(tsvL):
        # '.tags' middle term handling for groseq data vs normal data
        binName = path.basename(tsv).replace('.tags', '').replace('.tsv', ".bin")
        chrStr = binName.split('.')[0].replace('.bin', '')
        totalTagCnt = totalTagCnt + \
                      TSV2BIN(tsv, path.join(desDir,binName), \
                              maxClonalCnt, binsize, chrLenD[chrStr], fragLen, shiftMode)
    print totalTagCnt, maxClonalCnt
    fout = open(path.join(desDir, 'tagInfo.txt'), 'w')
    fout.write('TotalTagCount=%.1f\n' % totalTagCnt)
    fout.write('MaxClonalCount=%.1f\n' % maxClonalCnt)
    fout.close()
##    return totalTagCnt
        

# command line arguments:
# sys.argv[1]: srcDir or srcFile
# sys.argv[2]: desDir or srcFile
# sys.argv[3]: fragment length (optional)
if __name__=='__main__':
    
    genome, binsize, clonal_mode, fragLen = pipeUtil.readConfigSGE()
    chrLenD = genomeUtil.getChrLenD(genome)

    src, des = sys.argv[1], sys.argv[2]
    if len(sys.argv)>3: fragLen = sys.argv[3]
    if fragLen==0: fragLen = None   # fragment length not defined, use homer result
    if path.isdir(src):
        totalTagCnt = TSV2BIN_DIR(src, des, clonal_mode, binsize, chrLenD, fragLen)
        outDir = des
    else:
        pass #to be implemented
##        totalTagCnt = TSV2BIN(src, des, clonal_mode, binsize, chrLenD)
##        outDir = path.dirname(des)


    
