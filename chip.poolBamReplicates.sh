#!/usr/bin/env bash

###########################################3
# Cut&Run tools
# Written by Hee-Wooong Lim
# 
# Pool bam files of multiple replicates by group

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) <sample.tsv> <src bam directory> <des bam directory>
Description:
	Merge multiple bam files of a group according to sample/group information within a given sample.tsv file
	This script create a separate folder to place align.bam/log file per sample  unlike Cutlery verion, cnr.poolBamReplicates.sh
Input:
	- sample.tsv file: containing columns 'Name' and 'Group'
	- src bam directory: contanin sample folder containing align.bam file
	- des bam directory: to write sample folder containing a merged bam
Output:
	<des bam dir>/<sample name>/align.bam/log
Options:
	-b: if set, bsub are submitted for merging bam files, default=off
	-o: if set, overwrite existing destination bam files, default=off" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
bsub=FALSE
overwrite=FALSE
while getopts ":bo" opt; do
	case $opt in
		b)
			bsub=TRUE
			;;
		o)
			overwrite=TRUE
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
if [ $# -lt 3 ];then
	printUsage
	exit 1
fi
sampleInfo=$1
srcDir=$2
desDir=$3

assertFileExist $sampleInfo
assertDirExist $srcDir

###################################
## main code

#Usage: ngs.concateBamFiles.sh (options) [bam1] ...
#Description:
#	Concatenate bam files. No sorting
#Options:
#	-o <out>: Output file prefix including pathr. required
#	-m <mem>: Memory for sorting. default=5G
#	-s : If set, coordinate-sorted and indexed

groupL=`tail -n +2 $sampleInfo | grep -v -e ^$ -e ^# | cut -f 3 | sort | uniq`

echo -e "Pooling replicate bam files" >&2
echo -e "  - sampleInfo: $sampleInfo" >&2
echo -e "  - srcDir: $srcDir" >&2
echo -e "  - desDir: $desDir" >&2
echo -e "  - bsub:   $bsub" >&2

mkdir -p $desDir
for group in ${groupL[@]}
do
	des=${desDir}/${group}/align.bam
	log=${desDir}/${group}/align.log
	mkdir -p ${desDir}/${group}

	## Checking existing destination file
	if [ -f $des ] && [ "$overwrite" != "TRUE" ];then
		echo -e "Warning: $des already exists, pass" >&2
		continue
	fi

	## List of replicate bam files
	srcL=( `tail -n +2 $sampleInfo | grep -v -e ^$ -e "^#" | gawk '{ if($3 == "'$group'") printf "'$srcDir'/%s/align.bam\n", $2 }'` )
	assertFileExist ${srcL[@]}

	## Merging to destination
	echo -e "Creating $des" >&2
	for src in ${srcL[@]}
	do
		echo -e "  - $src" >&2
	done 2>&1 | tee $log

	if [ "$bsub" == "TRUE" ];then
		## Parallel processing using HPC:lsf
		bsub -W 24:00 -n 1 "module load samtools/1.9.0; ngs.concateBamFiles.sh -o $des ${srcL[@]}"
	else
		## Sequential processing
		ngs.concateBamFiles.sh -o $des ${srcL[@]}
	fi
done

