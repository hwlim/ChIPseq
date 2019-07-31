#/usr/bin/env bash

# Convert bin or nrm directory into a big wig file
# Usage:
#	BIN2BW.sh [srcDir] [desFile] [binsize]

if [ $# -lt 3 ];then
	echo -e "Usage: BIN2BW.sh [srcDir] [desFile] [binsize]"
	echo -e "TO BE IMPLEMENTED"
	exit
fi

src=$1
des=$2
binsize=$3

