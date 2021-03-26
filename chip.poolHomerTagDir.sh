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
	Merge multiple Homer tag directories of a group according to sample/group information within a given sample.tsv file
	Creating one tag directory per group
Input:
	- sample.tsv file: containing columns 'Name' and 'Group'
	- src sample directory: contanining a Homer Tag directory, TSV, i.e., <src sample dir>/<sampleName>/TSV
	- des group directory: destination sample directory to put sample folders, i.e., <des group dir>/<groupName>/TSV
Options:
	-b: if set, bsub are submitted for merging bam files, default=off
	-r : If set, perform robust estimation of fragment length, default=OFF
		because Homer gives unreliable fragment length when fragment length is close to the read length" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
bsub=FALSE
robust=FALSE
while getopts ":br" opt; do
	case $opt in
		b)
			bsub=TRUE
			;;
		r)
			robust=TRUE
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
if [ "$robust" = "TRUE" ];then
	optStr="-r"
else
	optStr=""
fi

groupL=`tail -n +2 $sampleInfo | grep -v -e ^$ -e ^# | cut -f 3 | sort | uniq`

echo -e "Pooling replicate bam files" >&2
echo -e "  - sampleInfo: $sampleInfo" >&2
echo -e "  - srcDir: $srcDir" >&2
echo -e "  - desDir: $desDir" >&2
echo -e "  - bsub:   $bsub" >&2
echo -e "  - robust: $robust" >&2
echo -e "" >&2

mkdir -p $desDir
for group in ${groupL[@]}
do
	des=${desDir}/${group}/TSV
	log=${desDir}/${group}/pool.log

	## List of replicate tag directories
	srcL=( `tail -n +2 $sampleInfo | grep -v -e ^$ -e "^#" | gawk '{ if($3 == "'$group'") printf "'$srcDir'/%s/TSV\n", $2 }'` )
	assertDirExist ${srcL[@]}

	## info.txt creation
	mkdir -p $des
	echo -e "Creating ${des}/../info.txt" >&2
	echo ${group} > ${des}/../info.txt

	## Merging to destination
	echo -e "Creating $des" >&2
	for src in ${srcL[@]}
	do
		echo -e "  - $src" >&2
	done 2>&1 | tee $log

	#continue
	## Pooling
	if [ "$bsub" == "TRUE" ];then
		## Parallel processing using HPC:lsf
		bsub -W 24:00 -M 1000 -n 1 "module load samtools/1.9.0; module load R/3.5.0; module load homer/4.11;
		ngs.poolHomerTagDir.sh -o $des $optStr ${srcL[@]};
		drawHomerAutoCorr.r -t ${group} ${des};"
	else
		## Sequential processing
		ngs.poolHomerTagDir.sh -o $des $optStr ${srcL[@]}
		drawHomerAutoCorr.r -t ${group} ${des}
	fi
done

