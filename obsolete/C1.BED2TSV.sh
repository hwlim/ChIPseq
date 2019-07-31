#!/usr/bin/env bash

export LC_ALL=C
LOGPATH="./log"
SRCDIR="./BED"
DESDIR="./TSV"

mkdir -p log
mkdir -p TSV

rm $LOGPATH/*.*

for src in $(ls $SRCDIR/*.bed)
do
	dataname=${src##*\/}
	dataname=${dataname%%.bed}
	HomerDes=$DESDIR/$dataname

	echo Submitting job:$src
#	echo HomerDestination:$HomerDes
	qsub -cwd -V -o $LOGPATH\/${dataname}.out -e $LOGPATH\/${dataname}.err ../scripts/BED2TSV.sh $src $HomerDes
done
unset LC_ALL
