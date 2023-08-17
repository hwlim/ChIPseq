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

grep ^Id $tsv
cat $tsv | grep -v -e ^$ -e "^#" -e ^Id \
	| gawk -F "\t" '{
				if(NF>7){
					printf "%s\t%s\t%s\n", $3,$7,$8
				}else{
					printf "%s\t%s\n", $3,$7
				}
			}' \
	| sort -k1,1 | uniq \
	| gawk -F "\t" '{
				if(NF>2){
					printf "%s\t%s\t%s\tNULL\tNULL\tNULL\t%s\t%s\n", $1, $1, $1, $2, $3
				}else{
					printf "%s\t%s\t%s\tNULL\tNULL\tNULL\t%s\n", $1, $1, $1, $2
				}
			}' 

