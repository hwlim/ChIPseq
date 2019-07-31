#!/usr/bin/env python

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
    

def gsChrBIN2WIG(srcFile, fout, binSize):
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
            


def gsBIN2WIG(srcFileL, desFile, trackName, binSize):
    fout = open(desFile, 'w')
    if desFile.find('.for.')!=-1: color = '255,0,0'
    else: color = '0,0,255'
    
    fout.write('track type=wiggle_0 name=%s description=%s autoScale=off viewLimits=0:100 visibility=full color=%s maxHeightPixels=100:50:30 windowingFunctoin=mean\n'\
               % (trackName, trackName, color))
    for srcFile in srcFileL:
##        print srcFile
        if srcFile.find('chrM')!=-1: continue
        print srcFile
        gsChrBIN2WIG(srcFile, fout, binSize)
    fout.close()

    

# arguments
#   srcDir, desDir, binSize
#   srcDir/*.for.bin(*.rev.bin) --> desDir/srcDir.for.wig(srcDir.rev.wig)
if __name__ == '__main__':
    if len(sys.argv)<4:
        print 'Arguments: srcDir, desDir, binSize'
        sys.exit(1)
    print sys.argv
    
    srcDir, desDir, binSize = \
            path.normpath(sys.argv[1]), path.normpath(sys.argv[2]), int(sys.argv[3])
    dataName = path.basename(srcDir)

    srcForL = glob.glob(path.join(srcDir, '*.for.bin'))
    srcRevL = glob.glob(path.join(srcDir, '*.rev.bin'))

    if len(srcForL) > 0 or len(srcRevL) > 0:
        desForFile = path.normpath(path.join(desDir, '%s.for.wig' % dataName))
        print desForFile
        gsBIN2WIG(sorted(srcForL), desForFile, dataName+'.forward', binSize)
        desRevFile = path.normpath(path.join(desDir, '%s.rev.wig' % dataName))
        print desRevFile
        gsBIN2WIG(sorted(srcRevL), desRevFile, dataName+'.reverse', binSize)
        os.system('gzip %s' % desForFile)
        os.system('gzip %s' % desRevFile)
    else:
        srcFileL = glob.glob(path.join(srcDir, '*.bin'))
        desFile = path.normpath(path.join(desDir, '%s.wig' % dataName))
        print desFile
        gsBIN2WIG(sorted(srcFileL), desFile, dataName, binSize)
        os.system('gzip %s' % desFile)

