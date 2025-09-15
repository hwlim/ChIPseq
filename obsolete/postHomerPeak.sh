#!/usr/bin/env bash

# Written by Hee Woong Lim, 12/15/2011
#
# Run homer with given options and perform post processing
# Usage:
#	run.homer.sh chipDir inputDir name
#	
#	chipDir/inputDir: chip and input tag directory
#	name: data name, all the results are generated with this name

source $COMMON_LIB_BASE/commonBash.sh

function printUsage {
	echo "Usage: postHomerPeak.sh [homerPeakFile] [tDir]" >&2
}

if [ $# -ne 2 ];then
	printUsage
	exit 1
fi

assertFileExist $1
isDirExist $2
peakFile=$1
tDir=$2
#genome=$3
outDir=`dirname $peakFile`
name=`basename $peakFile`
name=${name%%.txt}

grep -v '^#' $peakFile\
	| gawk -F'\t' '{printf "%s\t%s\t%s\t%s\t%.2f\t.\n", $2, $3, $4, $1, $6/10.0}'\
	> $outDir/${name}.bed

# 1rpm cut off, original
gawk '{if($5>=1){printf "%s\t%d\t%d\t%s\t%.2f\t%s\n", $1,$2,$3,$4,$5,$6}}' $outDir/${name}.bed \
	> $outDir/${name}.1rpm.bed

# 200bp width peak generation
gawk '{
	c=($2+$3)/2;
	if(c-100<=0)next;
	printf "%s\t%d\t%d\t%s\t%.2f\t%s\n", $1,c-100,c+100,$4,$5,$6;
	}' $outDir/${name}.bed \
	> ${outDir}/__temp__${name}.200bp.bed

tbp=""
ttc=""
tbp=`grep -m 1 "maximum tags considered per bp =" $peakFile | gawk '{print $(NF)}'`
ttc=`grep -m 1 "Tags Used for cluster" ${peakFile/%.txt/.log} | cut -d "=" -f 2 | gawk '{printf "%d", $1}'`
fLen=`grep -m 1 "fragment length" ${peakFile} | cut -d "=" -f 2 | gawk '{printf "%d", $1}'`
if [ "$tbp" == "" ];then
	echo "Error: Can't find tbp information" >&2
	exit 1
fi
if [ "$ttc" == "" ];then
	echo "Error: Can't find total tag count information" >&2
	exit 1
fi

#countTagsBed.sh -s both -m $tbp -na ${outDir}/__temp__${name}.200bp.bed $tDir\
countTagsBed.sh -s both -t $ttc -m $tbp -f $fLen -na ${outDir}/__temp__${name}.200bp.bed $tDir\
	| gawk '{printf "%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,$7,$6}'\
	| sort -k5,5nr\
	> ${outDir}/${name}.200bp.bed
rm ${outDir}/__temp__${name}.200bp.bed

# 1rpm cut off, 200bp
gawk '{if($5>=1){printf "%s\t%d\t%d\t%s\t%.2f\t%s\n", $1,$2,$3,$4,$5,$6}}' $outDir/${name}.200bp.bed \
	> ${outDir}/${name}.200bp.1rpm.bed


#echo "DONE!" >&2
#nPeak=`wc -l $outDir/${name}.200bp.1rpm.bed | cut -d " " -f 1`
#echo -e "Final number of peak calls(>=1rpm): $nPeak" >&2
