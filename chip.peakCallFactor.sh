#!/usr/bin/env bash

###########################################3
# Written by Hee-Wooong Lim
# 
# Wrapper script for peak calling
#
# - default option & additional option handling

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [taget tagDir]
Description: Make Homer data directory from BED file
Output:
	- <outDir>/peak.homer.txt              Homer peak calling result
	- <outDir>/peak.homer.bed              Homer peak in bed format
	- <outDir>/peak.homer.exBL.bed         After blacklist filtering
	- <outDir>/peak.homer.exBL.1rpm.bed    > 1rpm after filtering
Options:
	-o <outDir>: Destination tag directory, required
	-i <ctrl>: (optional) ctrl homer tag directory, default=NULL
	-m <mask>: mask bed file for filtering such as ENCODE blacklist
	-s <optStr>: additional option for 'findPeaks' of Homer
		Internal pre-set option: \"-style factor -tbp 0 -norm 1000000 -center\"
		Additional options are also possible such as -size 200 -minDist 400"
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
desDir=NULL
ctrl=NULL
mask=NULL
optStr=""
#lengthParam=NULL,NULL
while getopts ":o:i:m:s:" opt; do
	case $opt in
		o)
			desDir=$OPTARG
			;;
		i)
			ctrl=$OPTARG
			;;
		m)
			mask=$OPTARG
			;;
		s)
			optStr=$OPTARG
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

target=$1
assertDirExist $target

if [ "$desDir" == "NULL" ];then
	echo -e "Error: Destination directory (-o) must be specified" >&2
	exit 1
fi

if [ "$ctrl" != "NULL" ];then
	assertDirExist $ctrl
fi

if [ "$mask" != "NULL" ];then
	assertFileExist $mask
fi

###################################
## main code
log=${desDir}/peak.homer.log
echo -e "Homer peak-calling" >&2
echo -e "- target = $target" >&2
echo -e "- ctrl = $ctrl" >&2
echo -e "- desDir = $desDir" >&2
echo -e "- optStr = $optStr" >&2
echo -e "" >&2

peak0=${desDir}/peak.txt
peakBed=${desDir}/peak.bed
peakMasked=${desDir}/peak.exBL.bed
peak1rpm=${desDir}/peak.exBL.1rpm.bed

mkdir -p $desDir
if [ "$ctrl" == "NULL" ];then
	findPeaks $target -o ${peak0} -style factor -tbp 0 -norm 1000000 -center ${optStr} 2>&1 | tee ${log}
else
	findPeaks $target -i ${ctrl} -o ${peak0} -style factor -tbp 0 -norm 1000000 -center ${optStr} 2>&1 | tee ${log}
fi


grep -v "^#" ${peak0} \
	| gawk '{ printf "%s\t%d\t%d\t%s\t%s\t+\n", $2, $3, $4, $1, $6 }' \
	> ${peakBed}

if [ "$mask" != "NULL" ];then
	intersectBed -a ${peakBed} -b $mask -v > ${peakMasked}
else
	cp ${peakBed} ${peakMasked}
fi

gawk '$5 > 1' ${peakMasked} > ${peak1rpm}

N0=`cat ${peakBed} | wc -l`
N1=`cat ${peak1rpm} | wc -l`
echo -e "Final peaks: $N1 (/$N0), > 1rpm(/all)" >&2
