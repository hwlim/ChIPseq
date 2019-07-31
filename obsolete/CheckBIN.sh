#!/usr/bin/env bash

if [ $# -lt 1 ]; then
	echo -e "Usage:\n\tChecmkBIN.sh [BINDir]"
	exit
fi


export LC_ALL=C

BINPATH=$1
echo -e "DataName\tTotalTagCnt\tMaxClonalCount" > checkBIN.info
for bininfo in $(ls $BINPATH/*/tagInfo.txt)
do
	taginfo=`cat $bininfo`
	tagcnt=${taginfo##TotalTagCount=}
	maxclone=${tagcnt##*MaxClonalCount=}
	tagcnt=${tagcnt%%?MaxClonalCount=*}
	tagname=${bininfo%%\/tagInfo.txt}
	tagname=${tagname##*\/}
	echo -e "$tagname\t$tagcnt\t$maxclone" >> checkBIN.info
done
unset LC_ALL
