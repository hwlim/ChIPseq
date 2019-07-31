#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh

if [ $# -lt 2 ], then
	echo -e "Usage:\n\tBED2TSV.sh [tsvDIR] [bedFILE]"
	echo -e "ex) BED2TSV.sh TSV/dataname BED/bedfile.bed"
fi

isFileExist $2

makeTagDirectory $2 $1 -format bed

