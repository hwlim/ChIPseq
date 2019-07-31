#!/usr/bin/env python
#$ -S /gpfs/fs0/share/python2.6.5/Python-2.6.5/python
import sys, string, os, glob
import pipeUtil
import util
#from guppy import hpy
##import operator


##def Sub(H,I):
##    O = H[:]
##    for i in xrange(len(O)):
##        O[i] = H[i]-I[i]
##        #if O[i]<0: O[i]=0  # from Nha's suggestion
##    return O
##
##def Div(H,I):
##    O = H[:]
##    for i in xrange(len(H)):
##        O[i] = (H[i]+5.)/(I[i]+5.) -1.0
##        if O[i]<0: O[i]=0
##    return O

##def Average(Total):
##	denom = float(len(Total))
##	if denom==1: return Total[0]  
##	O = []
##	for i in xrange(len(Total[0])):
##		avg = sum(map(operator.itemgetter(i),Total))/denom
##		O.append(avg)
##	return O


##def PrintHistone(T,fname):
##    print "writing", fname
##    fp = open(fname,"w")
##    for i in xrange(len(T)):
##        fp.write("%.3f\n"%(T[i]))
##    fp.close()
    

def FGvsBG(fgNameL, bgName, srcDir, desDir):
#    h=hpy()
    for fgName in fgNameL:
        os.system('mkdir -p %s' % os.path.join(desDir, fgName))
#    h.heap()
        
    bgChrAll = glob.glob(os.path.join(srcDir, bgName, '*.nrm'))
    #print bgChrAll
    for bgChr in sorted(bgChrAll):
        bgNRM = util.readOneCol(bgChr, float)
#        h.heap()
        for fgName in fgNameL:
            fgChr = bgChr.replace(bgName, fgName)
            subPath = fgChr.replace(srcDir, desDir)
            print 'Writing %s' % subPath
            fin = open(fgChr, 'r')
            fout = open(subPath, 'w')
            for i, line in enumerate(fin):
                fout.write('%.3f\n' % (float(line)-bgNRM[i]))
            fin.close()
            fout.close()
##            data = util.readOneCol(fgChr, float)
###            h.heap()
##            subPath = fgChr.replace(srcDir, desDir)
##            data = Sub(data,bgNRM)
##            PrintHistone(data,subPath)
##            divName = subName.replace(".nrm", ".div")
##            DivH = Div(nrmFG,nrmBG)
##            PrintHistone(DivH,divName)



# sys.argv[1]: bgName
# sys.argv[2:]: a list of fgNames
if __name__=='__main__':
    bgName = sys.argv[1]
    fgNameL = []
    for i in xrange(len(sys.argv)-2): fgNameL.append(sys.argv[i+2])
    #print bgName
    #print fgNameL
    FGvsBG(fgNameL, bgName, './NRM', './SUB')

    

