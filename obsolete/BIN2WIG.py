#!/usr/bin/env python
# convert a give BIN/NRM directory into a wig file
# general version for GROseq and regular BIN/NRM files

import glob, os, string, sys
from os import path

class BINBuffer:
    def __init__(self, srcFile):
        self.fin = open(srcFile, 'r')
        self.current = float(self.fin.readline())
##        self.cnt = 1
        
    def __del__(self):
        self.fin.close()

    def readBlock(self):
        if self.current != None: block = [float(self.current)]
        else: return None

        while True:
            newLine = self.fin.readline()
            if newLine: self.current = float(newLine)
            else: self.current = None; break            
            if block[0] == 0:
                if self.current != 0: break
                else: block.append(self.current)
            else:
                if self.current == 0: break
                else: block.append(self.current)
        return block
    

def ChrBIN2WIG(srcFile, fout, binSize):
    zeroSpan = 50
    fileName = path.basename(srcFile)
    chrStr = fileName.split('.')[0]
##    if chrStr != 'chr2': return
    BINBuf = BINBuffer(srcFile)
    start = 1
    zeroState = False
    fmtStr = lambda x:'%d'%x
    
    while True:
        block = BINBuf.readBlock()
        if not block: break
        
        # fixtstep header output
        if start==1:
            if block[0]==0 and len(block)>=zeroSpan: # --> zero state
                zeroState = True
##                fout.write('fixedStep chrom=%s start=%d step=%d span=%d\n'\
##                           % (chrStr, start, len(block)*binSize, len(block)*binSize))
                fout.write('fixedStep chrom=%s start=%d step=%d\n'\
                           % (chrStr, start, len(block)*binSize))
                fout.write('0\n')
            else:
                fout.write('fixedStep chrom=%s start=%d step=%d\n'\
                           % (chrStr, start, binSize))
                fout.write('\n'.join(map(fmtStr, block))+'\n')
        else:
            if zeroState:
                zeroState=False
                fout.write('fixedStep chrom=%s start=%d step=%d\n'\
                           % (chrStr, start, binSize))
                fout.write('\n'.join(map(fmtStr, block))+'\n')
            else:
                if block[0]==0 and len(block)>=zeroSpan:
                    zeroState = True
##                    fout.write('fixedStep chrom=%s start=%d step=%d span=%d\n'\
##                               % (chrStr, start, len(block)*binSize, len(block)*binSize))
                    fout.write('fixedStep chrom=%s start=%d step=%d\n'\
                               % (chrStr, start, len(block)*binSize))
                    fout.write('0\n')
                else:
                    fout.write('\n'.join(map(fmtStr, block))+'\n')
##        fout.write('______________ zeroState=%s\n' % zeroState)
        start = start + binSize*len(block)
            


def BIN2WIG(srcDir, desDir, binSize, trackName, viewLimits, height):
    dataName = path.basename(srcDir)
    if not trackName: trackName = dataName
    if not viewLimits: viewLimits = 200
    if not height: height = 50
    
    srcFileL = glob.glob(path.join(srcDir, 'chr*.*'))
    desFile = path.normpath(path.join(desDir, '%s.wig' % dataName))
    
    fout = open(desFile, 'w')
    color = '0,0,255'    
    fout.write('track type=wiggle_0 name=%s description=%s autoScale=off viewLimits=0:%d visibility=full color=%s maxHeightPixels=100:%d:30 windowingFunctoin=mean\n'\
               % (trackName, trackName, viewLimits, color, height))
    for srcFile in srcFileL:
##        print srcFile
        if srcFile.find('chrM')!=-1: continue
##        print srcFile
        ChrBIN2WIG(srcFile, fout, binSize)
    fout.close()

    

# arguments
#   srcDir, desDir, binSize
#   srcDir: directory where a target bin files are located, data directory
#   desDir: directory where a wig file is generated, file name is the above srcDir
if __name__ == '__main__':
    if len(sys.argv)<4:
        print 'Arguments: srcDir, desDir, binSize, (trackName), (viewLimits), (height)'
        sys.exit(1)
    if len(sys.argv)<7: height=None
    else: height = int(sys.argv[6])
    if len(sys.argv)<6: viewLimits=None
    else: viewLimits = int(sys.argv[5])
    if len(sys.argv)<5: trackName=None
    else: trackName = int(sys.argv[4])

    ##    print sys.argv
    
    srcDir, desDir, binSize = \
            path.normpath(sys.argv[1]), path.normpath(sys.argv[2]), int(sys.argv[3])
    dataName = path.basename(srcDir)

    
    BIN2WIG(srcDir, desDir, binSize, trackName, viewLimits, height)
    os.system('gzip %s' % path.join(desDir, dataName+'.wig'))

