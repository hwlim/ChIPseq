#!/usr/bin/env bash

# command-line application template

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [bam]" >&2
	echo -e "Description:" >&2
	echo -e "\tPerform down-sampling for a given bam file and print randomly selected N record in BAM format to /dev/stdio" >&2
	echo -e "Options:" >&2
        echo -e "\t-n <N>: Number of record to select" >&2
        echo -e "\t-f : *not implemented yet* string of SAM flag option" >&2
        echo -e "\t-v : verbose mode" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
selectN=""
verbose=0
while getopts ":n:v" opt; do
	case $opt in
		n)
			selectN=$OPTARG
			;;
		v)
			verbose=1
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
if [ $# -eq 0 ];then
	printUsage
	exit 1
fi

src=$1
assertFileExist $src


###################################
## main code

downSample(){
	samtools view -H $1
	samtools view $1 \
		| shuf -n $selectN 
}

if [ $verbose -gt 0 ];then
	echo -e "Down sampling a BAM file" >&2
	echo -e "  src = $src" >&2
	echo -e "  N = $selectN" >&2
fi
downSample $src \
	| samtools view -b \
	| samtools sort -m 5G -T __temp__.$$ -

