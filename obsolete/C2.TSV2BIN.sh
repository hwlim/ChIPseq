#!/usr/bin/env bash

export LC_ALL=C # for case sensitive sorting 
LOGPATH="./log"
SRCDIR="./TSV"
DESBASEDIR="./BIN"

mkdir -p BIN
rm $LOGPATH/*.*

# src is each directory containing multiple chromosome
for src in $(ls -d $SRCDIR/*/)
do
	src=${src%\/} # elimination of last '/'
	srcname=${src##*\/}
	des=$DESBASEDIR/$srcname
	mkdir -p $des

	echo Submitting job:$src
	qsub -cwd -V -o $LOGPATH\/${srcname}.out -e $LOGPATH\/${srcname}.err ../scripts/TSV2BIN.py $src $des
done
unset LC_ALL
