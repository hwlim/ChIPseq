#!/usr/bin/env python
#$ -S /gpfs/fs0/share/python2.6.5/Python-2.6.5/python
# GRO-seq BED --> TSV
# ASSUMPTION: bed file is already sorted by chr, direction (optionaly start position)
import sys, os, string, glob
from os import path
##import genomeUtil

# TSV format
# chrStr position direction(0/1, +/-) cnt
# buffer class for read-in bed data block by chromosome/direction
class bedBuffer:
    def __init__(self, srcFile):
        self.fin = open(srcFile, 'r')        
        self.current = None
        self.isEmpty = False

    def __del__(self):
        self.fin.close()

    def readBED(self):
        if self.isEmpty: return None, None, None
        if not self.current:
            line = self.fin.readline()
            field = line.split()
            # chr, start, end, (seq), cnt, direction
            self.current = [field[0], field[1], field[2], field[4], field[5]]
            # 0:chr 1:start 2:end 3:cnt 4:direction

        chrStr = self.current[0]
        direc = self.current[4]
        if direc == '+': block = [[int(self.current[1]), float(self.current[3])]]
        elif direc == '-': block = [[int(self.current[2]), float(self.current[3])]]
        else: raise RuntimeError('Invalid direction: %s' % direc)
        
        while True:
            line = self.fin.readline()
            if line != '':
                field = line.split()
                if field[0] != self.current[0] or field[5] != self.current[4]:
                    self.current = [field[0], field[1], field[2], field[4], field[5]]
                    break
                else:
                    if direc=='+': block.append([int(field[1]), float(field[4])])
                    elif direc=='-': block.append([int(field[2]), float(field[4])])
                    else: raise RuntimeError('Invalid direction: %s' % direc)
            else:
                self.isEmpty = True
                break

        block.sort(key=lambda x:x[0])
        return chrStr, direc, block

    
class tagBuffer:
    def __init__(self, tagL):
        self.tagL = tagL
        self.current = None
        self.isEmpty = False

##    def __del__(self):

    def readTagCnt(self):
        if self.isEmpty or not self.tagL: return None, None
        if not self.current:
            self.current = 0
        pos = self.tagL[self.current][0]
        cnt = self.tagL[self.current][1]
        
        while True:
            self.current = self.current + 1
            if self.current < len(self.tagL):
                if pos != self.tagL[self.current][0]: break
                else: cnt = cnt + self.tagL[self.current][1]
            else:
                self.isEmpty = True
                break

        return pos, cnt


    
def gsBED2TSV(srcFile, desDir):
    print 'GRO-Seq BED2TSV: %s --> %s/' % (srcFile, desDir)
    os.system('mkdir -p %s' % desDir)

    cntD={} #[0]*maxDistrib # read count distribution, cntL[i]: numberb i+1 reads
    Nf, Nr = 0, 0   # total number of forward/reverse read count
    Npos = 0    # total number of read position

    bedBuf = bedBuffer(srcFile)
    while True:
        chrStr, direc, block = bedBuf.readBED()
        if not chrStr: break
##        print chrStr, direc, len(block)
        if direc=='+':
            tsvName = chrStr+'.for.tsv'
            direc=0
            Nf = Nf + sum(b[1] for b in block)
        elif direc=='-':
            tsvName = chrStr+'.rev.tsv'
            direc=1
            Nr = Nr + sum(b[1] for b in block)
        else:
            raise RuntimeError('Invalid Direction: %s' % direc)
        
        fout = open(path.join(desDir, tsvName), 'w')
        tagBuf = tagBuffer(block)
        while True:
            pos, cnt = tagBuf.readTagCnt()
            if not pos: break
            Npos = Npos + 1
            if cntD.has_key(int(cnt)): cntD[int(cnt)] = cntD[int(cnt)] + 1
            else: cntD[int(cnt)] = 1
            fout.write('%s\t%d\t%d\t%.3f\n' % (chrStr, pos, direc, cnt))
        fout.close()

    fout = open(path.join(desDir, 'tagInfo.txt'), 'w')
    fout.write('TotalReadCount=%d\n' % (Nf+Nr))
    fout.write('ForwardReadCount=%d\n' % Nf)
    fout.write('ReverseReadCount=%d\n' % Nr)
    fout.write('NumberOfGenomePosition=%d\n' % Npos)
    fout.write('AverageReadPerPosition=%.3f\n' % (float(Nf+Nr)/Npos))
    fout.write('\n===TagDistribution===\n')
    fout.write('ClonalCount\tReadCount\n')
    for i in sorted(cntD):
        fout.write('%d\t%d\n' % (i, cntD[i]))
    fout.close()    



def gsBED2TSV_Dir(srcDir, desDir):
    os.system('mkdir -p %s' % desDir)
    fileL = glob.glob(path.join(srcDir, '*.bed'))
    for bedFile in sorted(fileL):
        dataName = path.basename(bedFile).replace('.bed', '')
        gsBED2TSV(bedFile, path.join(desDir, dataName))
        



# command line arguments:
# sys.argv[1]: srcDir or srcFile
# sys.argv[2]: desDir
if __name__=='__main__':
    
    if len(sys.argv)<3:
        print 'Arguments: src, des'
        sys.exit(1)
    src, des = sys.argv[1], sys.argv[2]

    if path.isdir(src):
        gsBED2TSV_Dir(src, des)
    else:
        gsBED2TSV(src, des)
    
