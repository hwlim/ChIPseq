#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh
#trap 'rm __temp__.$$.*' EXIT

function printUsage {
	echo -e "Usage: picardMergeBam.sh <bamFile1> <bamFile2> ..." >&2
	echo -e "Description: Wrapper script for MarkDuplicates command in Picard tools" >&2
	echo -e "\tMerge & Sort give bam files and print-out to /dev/stdout" >&2
	echo -e "Options:" >&2
	echo -e "\t-s <sortOrder>: sort order, coordinate (default) / queryname / unsorted " >&2
	echo -e "\t-m <memory>: memory size, eg) 2000m (default), 2g ..." >&2

}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
sortOrder=coordinate
memory=2000m
while getopts ":s:m:" opt; do
	case $opt in
		s)
			sortOrder=$OPTARG
			;;
		m)
			memory=$OPTARG
			;;
		\?)
			echo "Invalid options: -$OPTARG" >&2
			printUsage
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			printUsage
			exit 1
			;;
	esac
done


shift $((OPTIND-1))
if [ $# -lt 2 ];then
	echo -e "Error: At least two bam files required for merging" >&2
	printUsage
	exit 1
fi
srcL=( $@ )
assertFileExist ${srcL[@]}



isInteger ${memory%[kmg]}


optStr=""
echo -e "Merging bam files by Picard" >&2
echo -e "  Src:" >&2
for src in ${srcL[@]}
do	
	echo -e "\t$src" >&2
	optStr="$optStr I=$src"
done
java -Xmx${memory} -jar ~/Programs/Picard/picard.jar MergeSamFiles \
	${optStr} \
	O=/dev/stdout \
	SORT_ORDER=${sortOrder} \
	USE_THREADING=true



