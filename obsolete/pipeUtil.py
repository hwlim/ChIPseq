#!/usr/bin/env python
import os, sys, glob, string

def normHistoneName(h):
    if h[2]=='k': h = h.replace('k', 'K')
    elif h=='Dnase': h='DNase'
    return h
    

def normHistoneNameRev(h):
    if h[2]=='K': h = h.replace('K', 'k')
    elif h=='DNase': h='Dnase'
    return h


def readConfigSGE(LEVEL=0):
    configD={}
    if LEVEL==0:    # basic level: genome, binsize, clonal_mode
        fin = open('config.txt', 'r')
        for line in fin:
            field = line[:-1].split('=')
            if len(field)!=2: break;
            configD[field[0]]=field[1]
        genome = configD['genome']
        binsize = int(configD['binsize'])
        clonal_mode = int(configD['clonal_mode'])
        fragLen = int(configD['fragmentLength'])
        fin.close()
        return genome, binsize, clonal_mode, fragLen

    elif LEVEL==1:  # full level: genome binsize, clonal_mode, fgD
        fin = open('config.txt', 'r')
        for line in fin:
            field = line[:-1].split('=')
            if len(field)!=2: break;
            configD[field[0]]=field[1]
        genome = configD['genome']
        binsize = int(configD['binsize'])
        clonal_mode = int(configD['clonal_mode'])
        fragLen = int(configD['fragmentLength'])
        
        fgD={}
        fgS=set()
        fgOnlyL=[]
        for line in fin:
            nameL = line.split()
            if len(nameL) == 0 or line[0]=='#': continue
            if len(nameL) == 2:
                fg, bg = nameL[0], nameL[1]
                if fg in fgS: RuntimeError('FG data appears more than once: %s' % fg)
                else: fgS.add(fg)
                if not fgD.has_key(bg): fgD[bg] = []
                fgD[bg].append(fg)
            elif len(nameL) == 1:
                fg = nameL[0]
                if fg in fgS: RuntimeError('FG data appears more than once: %s' % fg)
                else: fgS.add(fg)
                fgOnlyL.append(fg)
            else:
                raise RuntimeError('Invalid FG/BG pair, "%s"' % line[:-1])
        fin.close()
        return genome, binsize, clonal_mode, fragLen, fgD, fgOnlyL
    
    else:
        raise RuntimeError('Invalid LEVEL: %d', LEVEL)
    

def readConfig(LEVEL=0):
    configD={}
    if LEVEL==0:    # basic level: genome, binsize, clonal_mode, fragmentLength
        fin = open('config.txt', 'r')
        for line in fin:
            field = line[:-1].split('=')
            if len(field)!=2: break;
            configD[field[0]]=field[1]
        genome = configD['genome']
        binsize = int(configD['binsize'])
        clonal_mode = int(configD['clonal_mode'])
        fragLen = int(configD['fragmentLength'])
        
        fin.close()
        return genome, binsize, clonal_mode, fragLen

    elif LEVEL==1:  # full level: genome binsize, clonal_mode, fgD
        fin = open('config.txt', 'r')
        for line in fin:
            field = line.split();
            if len(field)==4: break     # break condition to be changed
            field = field[0].split('#')
            if not field[0]: continue
            field = field[0].split('=')
##            if len(field)!=2: break
            configD[field[0]]=field[1]
        genome = configD['genome']
        binsize = int(configD['binsize'])
        clonal_mode = int(configD['clonal_mode'])
        fragLen = int(configD['fragmentLength'])
        shiftMode = int(configD['shiftMode'])
        fgD={}
        fgS=set()
        fgOnlyL=[]
        fragLenD = {}
##        if fragLen == -1: fragLenD = {}
##        else: fragLenD = None
        namePairL = []
        
        for line in fin:
            field = line.split()
            if len(field) == 0 or line[0]=='#': continue
            assert(len(field) in [3,4])
            
            fg, bg, name = field[0], field[1], field[2]
            if bg=='None': bg = None
            
            if fg in fgS: RuntimeError('FG data appears more than once: %s' % fg)
            else: fgS.add(fg)
            if not fgD.has_key(bg): fgD[bg] = []
            fgD[bg].append(fg)

            if fragLen==-1: fragLenD[fg] = int(field[3])    # use individual value
            elif fragLen==0: fragLenD[fg] = None            # use Homer estimation
            else: fragLenD[fg] = fragLen                    # use fixed given value
            
            namePairL.append([fg,name])
            
        fin.close()
        return genome, binsize, clonal_mode, fragLen, shiftMode, fgD, namePairL, fragLenD
    
    else:
        raise RuntimeError('Invalid LEVEL: %d', LEVEL)                
    

def readConfig2():
    configD={}
    commonOptD={}
    fin = open('config.txt', 'r')
    for line in fin:
        field = line.split();
        if len(field)==4: break     # break condition to be changed
        field = field[0].split('#')
        if not field[0]: continue
        field = field[0].split('=')
##            if len(field)!=2: break
        configD[field[0]]=field[1]
    commonOptD['genome'] = configD['genome']
    commonOptD['binsize'] = int(configD['binsize'])
    commonOptD['clonal_mode'] = int(configD['clonal_mode'])
    commonOptD['fragLen'] = int(configD['fragmentLength'])
    commonOptD['shiftMode'] = int(configD['shiftMode'])
    dataOptD={}
    fgD={}
    fgS=set()
    fgOnlyL=[]
    fragLenD = {}
##        if fragLen == -1: fragLenD = {}
##        else: fragLenD = None
    namePairL = []

    serial=0;
    for line in fin:
        field = line.split()
        if len(field) == 0 or line[0]=='#': continue
        assert(len(field) in [3,4])
        
        fg, bg, name = field[0], field[1], field[2]
        if len(field)==4: fragLen=int(field[3])
        else: fragLen = 0

        if dataOptD.has_key(fg):
            RuntimeError('FG data appears more than once: %s' % fg)
        optD = {}
        optD['bg']=bg;
        optD['name']=name;
        optD['fragLen']=fragLen;
        optD['serial']=serial;
        serial = serial + 1;
        dataOptD[fg] = optD
        
    return commonOptD, dataOptD


def readHomerPeakConfig(File):
    configL=[]
    fin = open(File, 'r')
    for line in fin:
        line = line.splitlines()[0]
        optD={}
        field = line.split("\t");
        if len(field)==0 or field[0]== '' or field[0][0]=='#': continue
        optD['chip']=field[0]
        optD['input']=field[1]
##        print field
        if len(field)>2 and field[2] != '' : optD['name']=field[2]
        if len(field)>3 and field[3] != '' : optD['outDir']=field[3]
##        print '-%s-' % field[3]
        if len(field)>4 and field[4] != '' : optD['optionStr']=field[4]        
        configL.append(optD)
    return configL
    
##def readConfig(fname=None):
##    if not fname: fname = 'config.txt'
##    genome = None
##    binsize = None
##    maxproc = None
##    fgL = []
##    bgL = []
##    
##    fin = open(fname, 'r')
##    genome = fin.readline().split()[1]
##    binsize = int(fin.readline().split()[1])
##    maxproc = int(fin.readline().split()[1])
##    monoclonal = int(fin.readline().split()[1])
##    for line in fin:
##        if line[0]=='#': continue
##        nameL = line.split()
####        print nameL
##        if len(nameL)<1: continue
##        fgL.append(nameL[0])
##        if len(nameL)==2: bgL.append(nameL[1])
##    fin.close()
##
##    return genome, binsize, maxproc, monoclonal, fgL, bgL
                

# read tagInfo.txt generated by Homer
def readTSVInfo(tsvDir):
    fin = open(os.path.join(tsvDir,'tagInfo.txt'))
    fin.readline()  # header skip
    field = fin.readline().split()
    Position = int(field[1])
    TTC = float(field[2]) # total tag count
    fLen = int(fin.readline().split('=')[1])    # fragment length
    fin.readline()  # peak size estimate skip
    fin.readline()  # tags per BP skip
    avgTPP = float(fin.readline().split('=')[1]) # average tags per position
    fin.close()

    return TTC, fLen, avgTPP, Position


# read total tag cnt information in BIN directory
def readBINInfo(binDir):
    fin = open(os.path.join(binDir,'tagInfo.txt'))
    totalTagCnt = float(fin.readline().split('=')[1]) # total tag count
    fin.close()

    return totalTagCnt

    
## wgEncodeBroadHistoneGm12878ControlStdAlnRep1.bam --> Gm12878, Control, 1
def splitName(dataName, mode='ENCODE'):
    if mode=='ENCODE':
        uInd = []
        for i in xrange(len(dataName)):
            if dataName[i].isupper(): uInd.append(i)
    ##    print uInd
        if dataName[uInd[2]:uInd[3]]=='Dnase': 
            h = 'Dnase'
            repN = int(dataName[uInd[5]:].split('.')[0][3:])
        else: 
            h = dataName[uInd[4]:uInd[5]]
            repN = int(dataName[uInd[7]:].split('.')[0][3:])
        c = dataName[uInd[3]:uInd[4]]
        return c, h, repN
    
    elif mode=='EPIATLAS':
        field = dataName.split('.')
        if len(field)!=4:
            raise RuntimeError('Invalid data name format: %s' % dataName)
        c, h, serial = field[1:4]
        return c, h, serial
        
    else:
        raise RuntimeError('Invalid mode: %s' % mode)



####gawk '{OFS="\t"; if(strtonum($5)>20) print $1,$2,$3,$4,$5,$6}'
####  wgEncodeBroadHistoneGm12878ControlStdAlnRep1.bed | head
##def BAM2BED(srcPath, desPath, MQScrTh):
####    bamName = bamfile.split('/')[-1]
####    bedfile = bamfile.replace("bam", "bed")
##    
####    c, h, repN = splitName(bamfile)
####    print c, h, repN
##    
##    cmd = ("bamToBed -i %s | gawk " % srcPath \
##           + "'{OFS=\"\t\"; if(strtonum($5)>%d) print $1,$2,$3,$4,$5,$6}'" %MQScrTh \
##           + " > %s" % desPath)
##    print cmd
##    os.system(cmd)
