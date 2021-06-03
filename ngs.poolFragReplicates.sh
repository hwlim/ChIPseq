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
	Merge multiple fragment bed files of a group according to sample/group information within a given sample.tsv file
Input:
	- sample.tsv file: containing columns 'Name' and 'Group'
	- src fragment directory: contanin replicate bed files
	- des fragment directory: to write merged bed files

	Note:
	- frgament files are assumed to have a common suffix, *.frag.bed.gz
	- final fragment files are named as <group>.frag.bed.gz

Options:
	-b: if set, bsub are submitted for merging bam files, default=off
		Submission may be limited by LSF setting of total job submission limit
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

groupL=`tail -n +2 $sampleInfo | grep -v -e ^$ -e ^# | cut -f 3 | sort | uniq`

echo -e "Pooling replicate fragment bed files" >&2
echo -e "  - sampleInfo: $sampleInfo" >&2
echo -e "  - srcDir: $srcDir" >&2
echo -e "  - desDir: $desDir" >&2
echo -e "  - bsub:   $bsub" >&2

mkdir -p $desDir
for group in ${groupL[@]}
do
	des=${desDir}/${group}.frag.bed.gz
	log=${desDir}/${group}.log

	## Checking existing destination file
	if [ -f $des ] && [ "$overwrite" != "TRUE" ];then
		echo -e "Warning: $des already exists, pass" >&2
		continue
	fi

	## List of replicate bam files
	srcL=( `tail -n +2 $sampleInfo | grep -v -e ^$ -e "^#" | gawk '{ if($3 == "'$group'") printf "'$srcDir'/%s.frag.bed.gz\n", $2 }'` )
	assertFileExist ${srcL[@]}

	## Merging to destination
	echo -e "Creating $des" >&2
	for src in ${srcL[@]}
	do
		echo -e "  - $src" >&2
	done 2>&1 | tee $log

	## command to uncompress
	cmd="sort -m -k1,1 -k2,2n -k3,3n"
	for src in ${srcL[@]}
	do
		cmd="${cmd} <( zcat $src )"
	done
	echo $cmd >&2

	if [ "$bsub" == "TRUE" ];then
		## Parallel processing using HPC:lsf
		bsub -W 24:00 -n 1 "eval $cmd | gzip > $des"
	else
		## Sequential processing
		eval $cmd | gzip > $des
		#ngs.concateBamFiles.sh -o $des ${srcL[@]}
	fi
	echo "" >&2
done

