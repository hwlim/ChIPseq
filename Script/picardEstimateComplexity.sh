#!/usr/bin/env bash


source $COMMON_LIB_BASE/commonBash.sh
#trap 'rm __temp__.$$.*' EXIT

function printUsage {
	echo "Usage: picardSort.sh (options) [inputBam]" >&2
	echo "Description: Sort bam file by given criteria" >&2
	echo -e "Options:" >&2
#	echo -e "\t-s <sortOrder>: sort order, queryname / coordinate (default)" >&2
	echo -e "\t-m <memory>: memory size, eg) 2000m (default), 2g ..." >&2

}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
#sortOrder=coordinate
memory=2000m
while getopts ":m:" opt; do
	case $opt in
#		s)
#			sortOrder=$OPTARG
#			;;
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
if [ $# -lt 1 ];then
	printUsage
	exit 1
fi
assertFileExist $1

isInteger ${memory%[kmg]}

#if [ "$sortOrder" == "coordinate" ];then
#       suffix=csort.bam
#elif [ "$sortOrder" == "queryname" ];then
#       suffix=qsort.bam
#else
#	echo "Error: invalid sortOrder. must be coordinate or queryname" >&2
#	printUsage
#	exit 1
#fi


src=$1
des=`basename $src`
des=${des%.bam}.picardComplexity.txt


echo -e "Coordinate sorting by Picard" >&2
echo -e "  Src = $src" >&2
echo -e "  Des = $des" >&2
echo -e "  SortBy = $sortOrder" >&2
echo -e "  Memory = $memory" >&2
java -Xmx${memory} -jar ${PICARD_HOME}/picard.jar EstimateLibraryComplexity \
	I=${src} \
       	O=${des} \


