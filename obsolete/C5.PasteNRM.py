#!/usr/bin/env python
import sys, string, os
import pipeUtil

genome, binsize, clonal_mode, fgD = pipeUtil.readConfigSGE(LEVEL=1)

srcDir = './SUB'
pasteDir = './Paste'
os.system('mkdir -p Paste')
os.system('rm ./log/*.*')

##tmpDir = os.environ.get('TEMP')
##qsubShPath = os.path.join(tmpDir, '_sge_qsub.sh')
for i, bg in enumerate(sorted(fgD)):
    cell = pipeUtil.splitName(bg)[0]
    fgL = sorted(fgD[bg])

##    qsubShPath = os.path.join(tmpDir, '_sge_qsub.sh'+str(i))
##    fout = open(qsubShPath, 'w')
##    fout.write('#!/usr/bin/env bash\n')
##    fout.write('python ../scripts/PasteNRM.py %s %d %s %s\n' %
##               (genome, binsize, pasteDir, \
##                ' '.join([os.path.join(srcDir, fg) for fg in fgL])))
##    fout.close()

    print 'Submitting a job: %s' % cell
    cmd = 'qsub -cwd -V -o ./log/%s.out -e ./log/%s.err ../scripts/PasteNRM.py %s %d %s %s %s %s'\
          % (cell, cell, genome, binsize, srcDir, pasteDir, cell, ' '.join(fgL))
    os.system(cmd)
