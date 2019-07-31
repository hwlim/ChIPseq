#!/usr/bin/env bash

export LC_ALL=C # for case sensitive sorting 
LOGPATH="./log"
SRCBASEDIR="./BIN"
DESBASEDIR="./NRM"

mkdir -p NRM
rm $LOGPATH/*.*

# src is each directory containing multiple chromosome
for src in $(ls -d $SRCBASEDIR/*/)
do
	src=${src%\/} # elimination of last '/'
	srcname=${src##*\/}
	des=$DESBASEDIR/$srcname
	mkdir -p $des

	echo Submitting job:$src
	qsub -cwd -V -o $LOGPATH\/${srcname}.out -e $LOGPATH\/${srcname}.err ../scripts/BIN2NRM.py $src $des
done
unset LC_ALL
