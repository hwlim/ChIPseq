#!/usr/bin/env bash

###########################################3
# Written by Hee-Wooong Lim
# 
# Wrapper script for peak calling
#
# To do & consider:
# - default option & additional option handling

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [taget tagDir]
Description: Make Homer data directory from BED file
Options:
	-o <outPrefix>: output prefix including path, required
	-i <ctrl>: (optional) ctrl homer tag directory, default=NULL
	-m <mask>: mask bed file for filtering such as ENCODE blacklist
	-f <foldchange>: foldchange cut off for against control sample, default=4
	-k <spikein>: spike-in factor to multiply to the foldchange cut off (-f). (target spikein cnt)/(Input spikein cnt), default= 1
	-s <optStr>: additional option for 'findPeaks' of Homer
		Internal pre-set option: \"-style factor -norm 1000000 -center\"
		Warning: -tbp 0 is implicitly set (may be revised in the future)
		Additional options are also possible such as -size 200 -minDist 400
Output:
	- <outPrefix>.txt              Homer peak calling result
	- <outPrefix>.bed              Homer peak in bed format
	- <outPrefix>.exBL.bed         After blacklist filtering
	- <outPrefix>.exBL.1rpm.bed    > 1rpm after filtering"
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
outPrefix=NULL
ctrl=NULL
mask=NULL
foldchange=4
spikein=1
optStr=""
#lengthParam=NULL,NULL
while getopts ":o:i:m:f:k:s:" opt; do
	case $opt in
		o)
			outPrefix=$OPTARG
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

if [ "$outPrefix" == "NULL" ];then
	echo -e "Error: outPrefix (-o) must be specified" >&2
	exit 1
fi

if [ "$ctrl" != "NULL" ];then
	assertDirExist $ctrl
fi

if [ "$mask" != "NULL" ];then
	assertFileExist $mask
fi

desDir=`dirname $outPrefix`
mkdir -p $desDir

###################################
## main code
echo -e "Homer TF peak-calling
  - target = $target
  - ctrl = $ctrl
  - blacklist = $mask
  - outPrefix = $outPrefix
  - fold-change = $foldchange
  - spikein = $spikein" >&2

log=${outPrefix}.homer.log
peak0=${outPrefix}.txt
peakBed=${outPrefix}.bed
peakMasked=${outPrefix}.exBL.bed
peak1rpm=${outPrefix}.exBL.1rpm.bed


## optional fold-change cut off adjustment by spikein factor
if [ "$spikein" != "1" ];then
	foldchange=`echo -e "${foldchange}\t${spikein}" | gawk '{ printf "%f", $1 * $2 }'`
	echo -e "  - Spikein factor applied, new fold-change threshold:
	spikein factor = $spikein
	new fold change = $foldchange
" >&2
fi
optStr="$optStr -F $foldchange"
echo -e "  - optStr = $optStr\n" >&2


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
