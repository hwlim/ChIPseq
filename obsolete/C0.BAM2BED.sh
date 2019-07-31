#!/usr/bin/env bash

export LC_ALL=C
mkdir -p log
LOGPATH="./log"
DATAPATH="./BAM"
DESDIR="./BED"

mkdir -p BED
rm $LOGPATH/*.*

for src in $(ls $DATAPATH/*.bam)
do
	dataname=${src##*\/}
	dataname=${dataname%%.bam}
	des=$DESDIR/${dataname}.bed

	echo "Submitting job:$src"
	qsub -cwd -V -o $LOGPATH\/${dataname}.out -e $LOGPATH\/${dataname}.err ../scripts/BAM2BED.sh $src $des $1
done
unset LC_ALL
