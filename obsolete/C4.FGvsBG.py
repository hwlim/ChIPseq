#!/usr/bin/env python
import sys, string, os
import pipeUtil

genome, binsize, clonal_mode, fragLen, fgbgD , fgOnlyL = pipeUtil.readConfigSGE(LEVEL=1)

os.system('mkdir -p SUB')
os.system('rm ./log/*.*')


for i, bg in enumerate(sorted(fgbgD)):
    cmd = 'qsub -cwd -V -o ./log/%s.out -e ./log/%s.err ../scripts/FGvsBG.py %s' \
          % (bg, bg, ' '.join([bg]+sorted(fgbgD[bg])))
##    print cmd
    print 'Submitting a job for a background: %s' % bg
    os.system(cmd)

for fg in fgOnlyL:
    print 'Copying fgOnly data to SUB: %s' % fg
    os.system('cp -r NRM/%s SUB/' % fg)
