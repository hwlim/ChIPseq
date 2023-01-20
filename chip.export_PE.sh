#!/usr/bin/env bash

###########################################3
# Written by Hee-Wooong Lim
# 

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) <src dir> <des dir>
Description: Export ChIP-seq/PE processing results to a destination folder
Input:
	A folder containing sample folders that containing outputs, i.e. stratified by samples
	e.g.
	3.Sample/
	├── <sample1>
	├── <sample2>
	└── <sample3>
Output:
	A folder containing outputs organized by output types and prefixed by sample name
	e.g.
	BigWig/
	├── <sample1>.bw
	├── <sample2>.bw
	└── <sample3>.bw
	Peak/
	├── <sample1>.bed
	├── <sample2>.bed
	└── <sample3>.bed

Options:
	-l : Create a symbolic link instead of copying for large files including
		bigwig and motif search result
	-f : Overwrite existing files/folders if set" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi



###################################
## option and input file handling
force=FALSE
makelink=FALSE
#lengthParam=NULL,NULL
while getopts ":lf" opt; do
	case $opt in
		l)
			makelink=TRUE
			;;
		f)
			force=TRUE
			;;
#		s)
#			optStr=$OPTARG
#			;;
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
	printUsage
	exit 1
fi

srcDir=$1
desDir=$2

assertDirExist $srcDir
mkdir -p $desDir


echo -e "Exporting Cutlery results" >&2
echo -e "  src: $srcDir" >&2
echo -e "  des: $desDir" >&2
echo -e "  force: $force" >&2
echo -e "  makelink: $makelink" >&2


sampleL=( `ls -d ${srcDir}/*/` )

## export list
# - bigwig: all / nfr / nuc
# - peak & heatmap: all / 1rpm
# - motif

exportFile(){
	local src=$1
	local des=$2
	local isDir=$3
	local makelink=$4

	if [ "$isDir" == "TRUE" ];then
		optStr="-r"
	else
		optStr=""
	fi

	if [ -f $src ] || [ -d $src ];then
		local desDir=`dirname $des`
		mkdir -p $desDir

		if [ -f $des ];then
			if [ "$force" == "TRUE" ];then
				echo -e "  Warning: Overwriting $des" >&2
				if [ "$makelink" == "TRUE" ];then
					rm -f $des
					ln -svr $src $des
				else
					cp $optStr -f -v $src $des
				fi
			else
				echo -e "  Warning: $des already exists, skip" >&2
			fi
		else
			if [ "$makelink" == "TRUE" ];then
				ln -svr $src $des
			else
				cp $optStr -v $src $des
			fi
		fi
	else
		echo -e "  Warning: $src does not exist, skipping" >&2
	fi
}

## Export selection
read -p 'Export QC result [Y/n]:' exportQC
read -p 'Export bigwig files [Y/n]:' exportBW
read -p 'Export peak files [Y/n]:' exportPeak
read -p 'Export motif results [Y/n]:' exportMotif
exportQC=${exportQC:-Y}
exportBW=${exportBW:-Y}
exportPeak=${exportPeak:-Y}
exportMotif=${exportMotif:-Y}

## Export
echo -e "" >&2
for samplePath in ${sampleL[@]}
do
	echo -e "Exporting ${samplePath}" >&2
	sample=`basename $samplePath`

	# fragment length distribution
	if [ "$exportQC" == "y" ] || [ "$exportQC" == "Y" ];then
		echo -e "0) Exporting fragment length distribution" >&2
		#mkdir -p ${desDir}/QualityControl
		src=${srcDir}/${sample}/QC/fragLen.png
		des=${desDir}/QualityControl/${sample}.fragLen.png
		exportFile $src $des FALSE FALSE
	else
		echo -e "0) Skipping fragment length distribution" >&2
	fi

	# - bigwig: all / nfr / nuc
	if [ "$exportBW" == "y" ] || [ "$exportBW" == "Y" ];then
		echo -e "1) Exporting bigwig files" >&2
		#mkdir -p ${desDir}/BigWig
		for prefix in "" .ctr.rpm .ctr.rpsm .frag.rpm .frag.rpsm
		do
			for suffix in "" ".subInput"
			do
				src=${srcDir}/${sample}/igv${prefix}${suffix}.bw
				des=${desDir}/BigWig${prefix}${suffix}/${sample}${prefix}${suffix}.bw
				exportFile $src $des FALSE $makelink
			done
		done
	else
		echo -e "1) Skipping bigwig files" >&2
	fi

	# - peak & heatmap: all / 1rpm
	if [ "$exportPeak" == "y" ] || [ "$exportPeak" == "Y" ];then
		echo -e "2) Exporting peak files & heatmaps" >&2
		#mkdir -p ${desDir}/Peak
		for mode in factor histone
		do
			src=${srcDir}/${sample}/HomerPeak.${mode}/peak.exBL.bed
			des=${desDir}/Peak.${mode}/${sample}.all.bed
			exportFile $src $des FALSE FALSE

			src=${srcDir}/${sample}/HomerPeak.${mode}/peak.exBL.1rpm.bed
			des=${desDir}/Peak.${mode}/${sample}.1rpm.bed
			exportFile $src $des FALSE FALSE

			src=${srcDir}/${sample}/HomerPeak.${mode}/heatmap.exBL.1rpm.png
			des=${desDir}/Peak.${mode}/${sample}.1rpm.heatmap.png
			exportFile $src $des FALSE FALSE
		done
	else
		echo -e "2) Skipping peak files & heatmaps" >&2
	fi

	# - motif
	if [ "$exportMotif" == "y" ] || [ "$exportMotif" == "Y" ];then
		echo -e "3) Exporting motif search results" >&2
		#mkdir -p ${desDir}/Motif
		for src in ${srcDir}/${sample}/HomerPeak.factor/Motif/Homer.all ${srcDir}/${sample}/HomerPeak.factor/Motif/MEME.random5k
		do
			suffix=`basename $src`
			des=${desDir}/Peak.factor/Motif/${sample}.${suffix}
			exportFile $src $des TRUE $makelink
		done
	else
		echo -e "3) Skipping motif search results" >&2
	fi
done
