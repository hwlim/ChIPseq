#!/usr/bin/env bash

#########################################
# Lim, Hee-Woong
# Jan/09/2017
#########################################


source $COMMON_LIB_BASE/commonBash.sh

if [ $# -lt 1 ];then
	echo -e "starClearMemory.sh <star index directory>" >&2
	echo -e "Description: Remove STAR index from a shared memory if it is not being used" >&2
	exit 0
fi

starIndex=$1
if [ ! -d $starIndex ];then
	echo -e "Error: $starIndex is not a valid directory" >&2
	exit 1
else
	starIndex=`readlink -m $starIndex`
fi


set +o pipefail
indexBeingUsed=`ps aux | grep -w STAR | grep $starIndex | grep -v grep | wc -l`
if [ $indexBeingUsed -gt 0 ];then
        echo -e "STAR index ${starIndex} is being used by other job, will be kept in a shared memory" >&2
else
        echo -e "STAR index ${starIndex} is not being used, will be removed from a shared memory" >&2
        STAR --genomeDir ${starIndex} --genomeLoad Remove
fi


