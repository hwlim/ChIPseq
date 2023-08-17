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
	-m <mask>: mask bed file for filtering such as ENCODE blacklist, default=NULL (no filtering)
	-s <optStr>: additional option for 'findPeaks' of Homer
		in addition to internally pre-set option: \"-style histone -tbp 0 -norm 1000000 \"
		such as -size or -minDist
		Warning: -tbp 0 is implicitly set (may be revised in the future)
Output:
	- <outPrefix>.txt              Homer peak calling result
	- <outPrefix>.bed              Homer peak in bed format
	- <outPrefix>.exBL.bed         After blacklist filtering
	- <outPrefix>.param.txt        peak calling parameter specified by -s"
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
optStr="-style histone -tbp 0 -norm 1000000 -strand both"
while getopts ":o:i:m:s:" opt; do
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
			optStr="$optStr $OPTARG"
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

assertFileExist ${target}/tagInfo.txt
ttc=`( grep genome ${target}/tagInfo.txt || true ) | cut -f 3`
if [ "$ttc" == "" ];then
	echo -e "Error: no tag count information found in [${target}/tagInfo.txt]" >&2
	exit 1
fi

###################################
## main code
echo -e "Homer histone peak-calling
  - target = $target
  - ctrl = $ctrl
  - blacklist = $mask
  - outPrefix = $outPrefix
  - TTC = $ttc
  - option = $optStr" >&2

log=${outPrefix}.log
param=${outPrefix}.param.txt
peak0=${outPrefix}.txt
peakBed=${outPrefix}.bed
peakMasked=${outPrefix}.exBL.bed

tmpPeakMasked=${TMPDIR}/__temp__.$$.bed
tmpTagCount=${TMPDIR}/__temp__.$$.target

echo -e "$optStr" > $param

desDir=`dirname $outPrefix`
mkdir -p $desDir
if [ "$ctrl" == "NULL" ];then
	echo -e "findPeaks $target -o ${peak0} ${optStr}" >&2
	findPeaks $target -o ${peak0} ${optStr} 2>&1 | tee ${log}
else
	echo -e "findPeaks $target -i ${ctrl} -o ${peak0} ${optStr}" >&2
	findPeaks $target -i ${ctrl} -o ${peak0} ${optStr} 2>&1 | tee ${log}
fi



echo -e "Convergting to bed" >&2
grep -v "^#" ${peak0} \
	| gawk '{ printf "%s\t%d\t%d\t%s\t%.3f\t+\n", $2, $3, $4, $1, $6*1000/($4-$3) }' \
	> ${peakBed}

echo -e "Blacking masking & merging" >&2
if [ "$mask" != "NULL" ];then
	#intersectBed -a ${peakBed} -b $mask -v \
	subtractBed -a ${peakBed} -b $mask  \
		| sortBed \
		| mergeBed \
		| gawk '$3-$2>=500 { printf "%s\t%d\t%d\tpeak.%d\t0\t+\n", $1,$2,$3,NR }' \
		> ${tmpPeakMasked}
		#| sort -k4,4 \
else
	cat ${peakBed} \
		| sortBed \
		| mergeBed \
		| gawk '{ printf "%s\t%d\t%d\tpeak.%d\t0\t+\n", $1,$2,$3,NR }' \
		> ${tmpPeakMasked}
fi

if [ `cat $tmpPeakMasked | wc -l` -eq 0 ];then
	echo -e "No remaining peaks" >&2
	touch $peakMasked
else
	echo -e "Tag counts in RPKM" >&2

	getPeakTags $tmpPeakMasked $target -tagAdjust 0 -tbp 0 -fixed \
		> ${tmpTagCount}

	# coordinate / line number / RPKM / strand / size / RPM
	paste <( sort -k4,4 ${tmpPeakMasked} ) <( sort -k1,1 ${tmpTagCount} ) \
		| gawk '{ printf "%s\t%d\t%d\tpeak.%d\t%.5f\t%s\t%d\t%.1f\n", $1,$2,$3, NR, $8*1000000/'${ttc}'*1000/($3-$2), $6,$3-$2, $8*1000000/'${ttc}' }' \
		| sort -k8,8nr \
		| cut -f 1-6 \
		> $peakMasked
fi

N0=`cat ${peakBed} | wc -l`
N1=`cat ${peakMasked} | wc -l`
echo -e "Final peaks:" >&2
echo -e "  - Original = $N0" >&2
echo -e "  - After filtering = $N1" >&2

