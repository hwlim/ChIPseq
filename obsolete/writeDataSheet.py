#!/usr/bin/env python
import os, sys, glob, string


srcPath = './BAM/'

def splitName(bamName):
    uInd = []
    for i in xrange(len(bamName)):
        if bamName[i].isupper(): uInd.append(i)
##    print uInd
    c = bamName[uInd[3]:uInd[4]]
    h = bamName[uInd[4]:uInd[5]]
    repN = int(bamName[uInd[7]:].split('.')[0][3:])

    return c, h, repN



bamDataL = glob.glob(srcPath+'*.bam')
bamDataL = [bamData.replace('\\', '/') for bamData in bamDataL] # window compatibility
bamDataL = [bamData.split('/')[-1].split('.')[0] for bamData in bamDataL]


hD = {}
cD = {}
hL = set()
cL = set()
for bamData in bamDataL:
    c, h,  repN = splitName(bamData)
    hL.add(h)
    cL.add(c)
    if not hD.has_key(c): hD[c] = set()
    hD[c].add(h)
    if not cD.has_key(h): cD[h] = set()
    cD[h].add(c)
hL = sorted(list(hL))
cL = sorted(list(cL))


os.system('mkdir -p Info')

fout = open('./Info/histoneInCell.txt', 'w')
for c in hD: 
    fout.write('\t'.join([c+':']+sorted(list(hD[c])))+'\n')
fout.close()

fout = open('./Info/cellInHistone.txt', 'w')
for h in cD:
    fout.write('\t'.join([h+':']+sorted(list(cD[h])))+'\n')
fout.close()



fout = open('./Info/cellHistoneTable.txt', 'w')
fout.write('Cells\t'+'\t'.join(hL)+'\n')
for c in cL: 
    fout.write('\t'.join([c+':']+['T' if h in hD[c] else 'F' for h in hL])+'\n')
fout.write('\n\n')
fout.write('Histons\t'+'\t'.join(cL)+'\n')
for h in hL: 
    fout.write('\t'.join([h+':']+['T' if c in cD[h] else 'F' for c in cL])+'\n')
fout.close()

for c in hD:
    fout = open('./Info/%s.histone'%c, 'w')
    fout.write('\n'.join(hD[c])+'\n')
    fout.close()

for h in cD:
    fout = open('./Info/%s.cell'%h, 'w')
    fout.write('\n'.join(cD[h])+'\n')
    fout.close()
