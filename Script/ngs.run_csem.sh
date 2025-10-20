#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.run_csem.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.run_csem.$$.*; fi' EXIT

function printUsage {
	echo "Usage: `basename $0` (options) [bam]
Description:
	Wrapper script to run CSEM
	This script was created to encapsulate average fragment length determination for Snakemake
Input:
	Name-sorted paired-end bam file containing multi-mappers
Options:
	-o : output prefix. required
	-t : number of thread
	-n : name for plot title, If not specified, bam file name is used" >&2
}

###################################
## option and input file handling
outPrefix=NULL
thread=1
name=NULL
while getopts ":o:t:n:" opt; do
	case $opt in
		o)
			outPrefix=$OPTARG
			;;
		t)
			thread=$OPTARG
			;;
		n)
			name=$OPTARG
			;;
		# c)
		# 	chrRegex=$OPTARG
		# 	;;
		# r)
		# 	robust=TRUE
		# 	;;
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

if [ "$outPrefix" = "NULL" ];then
	echo -e "Error: outPrefix (-o) must be specified" >&2
	exit 1
fi

if [ "$name" = "NULL" ];then
	name=$src
fi

echo -e "============================================
Running CSEM
============================================
  - src: $src
  - outPrefix: $outPrefix
  - name: $name
  - thread: $thread" >&2

desTxt=${outPrefix}.fragLen.txt

outDir=`dirname $outPrefix`
mkdir -p $outDir


echo -e "Calculating average fragment length" >&2
echo -e "Running
> ngs.fragLenHist.r -o ${outPrefix}.fragLen -l 1000 -n $name -a $src">&2
avg=`ngs.fragLenHist.r -o ${outPrefix}.fragLen -l 1000 -n $name -a $src`
if [ $avg -lt 10 ];then
	echo -e "Error: average lenght is shorter than 10 bp; something's wrong" >&2
	exit 1
fi
# ngs.getFragLenHist.sh $src > $desTxt
# avg=`ngs.getAvgFragLen.r ${desTxt}`
echo -e "  - ${avg} bp" >&2
echo -e "Running CSEM" >&2
#run-csem --bam -p ${thread} ${src} ${avg} ${outPrefix}
# csem b test.chr10.bam 150 200 test2 4 --extend-reads"
csem b ${src} ${avg} 200 ${outPrefix} ${thread}
