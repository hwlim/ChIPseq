#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh

if [ $# -lt 1 ];then
	echo -e "Usage: chip.make_pool_tsv.sh <sample.tsv>
Description:
	Create sample.tsv for pooled data processing" >&2
	exit 1
fi

tsv=$1

assertFileExist $tsv

head -n 1 $tsv
tail -n +2 $tsv | cut -f 3,7 | sort -k1,1 | uniq |  gawk '{ printf "%s\t%s\t%s\tNULL\tNULL\tNULL\t%s\n", $1, $1, $1, $2 }' 

