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
Input:
	- sample.tsv file: containing columns 'Name' and 'Group', from original/replicate Cutlery run
	- src bam directory: contaning sample folders containing replicate bam file
		Case 1: all bam files in a single foldler
			<src sample dir>
			├── <sample1>.bam
			└── <sample2>.bam
		Case 2: stratified by sample folder
			<src sample dir>
			├── sample1
			│    └── align.bam
			└── sample12
			     └── align.bam
	- des bam directory: to write merged bam files
		Case 1:
			<des group dir>
			├── <group1>.bam
			└── <group2>.bam
		Case 2:
			<des group dir>
			├── group1
			│    └── align.bam
			└── group2
			     └── align.bam
Options:
	-t: if set, print pooling plan lisintg src and des files
		No actual pooling job, just for simulation
	-s: if set, input bam files are stratified by sample folders as shown in Case 2
		In default, all input bam files are assumed to be in a single src folder as shown in Case 1
	-b: if set, bsub are submitted for merging bam files, default=off
	-u: if set, assume unsorted input bam files and create unsorted output, default=off
		NOTE: Previousu Cutlery created unsorted bam file. But the current version create coordinate-sorted bam files after alignment
	-f: if set, force overwrite existing destination bam files, default=off" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
bsub=FALSE
unsortedBam=FALSE
overwrite=FALSE
testOnly=FALSE
bySampleDir=FALSE
while getopts ":bsuft" opt; do
	case $opt in
		b)
			bsub=TRUE
			;;
		s)
			bySampleDir=TRUE
			;;
		u)
			unsortedBam=TRUE
			;;
		f)
			overwrite=TRUE
			;;
		t)
			testOnly=TRUE
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
echo -e "  - bySampleDir: $bySampleDir" >&2
echo -e "  - testOnly: $testOnly" >&2
echo -e "  - unsorted bam: $unsortedBam" >&2
echo -e "  - bsub:   $bsub" >&2
echo -e "" >&2

mkdir -p $desDir
for group in ${groupL[@]}
do
	if [ "$bySampleDir" == "TRUE" ];then
		des=${desDir}/${group}/align.bam
		log=${desDir}/${group}/align.log
	else
		des=${desDir}/${group}.bam
		log=${desDir}/${group}.log
	fi

	## Checking existing destination file
	if [ -f $des ] && [ "$overwrite" != "TRUE" ];then
		echo -e "Warning: $des already exists, pass" >&2
		continue
	fi

	## List of replicate bam files
	if [ "$bySampleDir" == "TRUE" ];then
		srcL=( `tail -n +2 $sampleInfo | grep -v -e ^$ -e "^#" | gawk '{ if($3 == "'$group'") printf "'$srcDir'/%s/align.bam\n", $2 }'` )
	else
		srcL=( `tail -n +2 $sampleInfo | grep -v -e ^$ -e "^#" | gawk '{ if($3 == "'$group'") printf "'$srcDir'/%s.bam\n", $2 }'` )
	fi
	assertFileExist ${srcL[@]}

	## Merging to destination
	echo -e "Creating $des" >&2
	for src in ${srcL[@]}
	do
		echo -e "  - $src" >&2
	done

	if [ "$testOnly" == "TRUE" ];then
		continue
	fi

	if [ "$unsortedBam" == "TRUE" ];then
		optStr=""
	else
		optStr="-S"
	fi

	if [ "$bySampleDir" == "TRUE" ];then
		mkdir -p ${desDir}/${group}
	fi

	echo -ne "" > $log
	for src in ${srcL[@]}
	do
		echo -e "- $src" >> $log
	done

	if [ "$bsub" == "TRUE" ];then
		## Parallel processing using HPC:lsf
		bsub -W 24:00 -n 1 "module load samtools/1.9.0; ngs.concateBamFiles.sh $optStr -o $des ${srcL[@]}"
	else
		## Sequential processing
		ngs.concateBamFiles.sh $optStr -o $des ${srcL[@]}
	fi
done

