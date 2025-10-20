#!/usr/bin/env bash

# command-line application template

source $COMMON_LIB_BASE/commonBash.sh
#trap 'rm __temp__.$$.*' EXIT

function printUsage {
	echo -e "Usage: getBowtieStats.sh (options) [inputFiles] ..." >&2
	echo -e "InputFile: bowtie log file" >&2
	echo -e "Options: <no options yet>" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
#while getopts ":g:f0" opt; do
#	case $opt in
#		g)
#			genome=$OPTARG
#			;;
#		f)
#			isFlag=1
#			;;
#		0)
#			echo -e "Flat -0 set" >&2
#			;;
#		\?)
#			echo "Invalid options: -$OPTARG" >&2
#			printUsage
#			exit 1
#			;;
#		:)
#			echo "Option -$OPTARG requires an argument." >&2
#			printUsage
#			exit 1
#			;;
#	esac
#done


shift $((OPTIND-1))
if [ $# -eq 0 ];then
	printUsage
	exit 1
fi

srcL=( $@ )
assertFileExist $@

###################################
## main code

echo -e "Data\tTotalReads\tAligned\tFailed\tSuppressed\tAligned%\tFailed%\tSuppressed%" 
for src in ${srcL[@]}
do
	name=`basename $src`
	name=${name%.log}

	ttc=`grep  "# reads processed" $src | cut -d " " -f 4`
	aligned=`grep "# reads with at least one" $src | cut -d " " -f 9,10`
	failed=`grep "# reads that failed" $src | cut -d " " -f 7,8`
	suppressed=`grep "# reads with alignments suppressed" $src | cut -d " " -f 9,10`

	align_N=`echo $aligned | cut -d " " -f 1`
	align_R=`echo $aligned | cut -d " " -f 2 | sed 's/[()]//g'`
	failed_N=`echo $failed | cut -d " " -f 1`
	failed_R=`echo $failed | cut -d " " -f 2 | sed 's/[()]//g'`
	sprsd_N=`echo $suppressed | cut -d " " -f 1`
	sprsd_R=`echo $suppressed | cut -d " " -f 2 | sed 's/[()]//g'`

	echo -e "${name}\t${ttc}\t${align_N}\t${failed_N}\t${sprsd_N}\t${align_R}\t${failed_R}\t${sprsd_R}"
done
