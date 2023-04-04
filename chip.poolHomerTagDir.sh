#!/usr/bin/env bash

###########################################3
# Cut&Run tools
# Written by Hee-Wooong Lim
# 
# Pool multiple homer tag dirs replicates by group

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) <sample.tsv> <src sample directory> <des sample directory>
Description:
	Merge multiple Homer tag directories of a group according to sample/group information within a given sample.tsv file
	Creat one tag directory per group
	** Skip already existing destination folder
	** If there is only one replicate, simply symbolic link is created not performing pooling.
Input:
	- sample.tsv file: containing columns 'Name' and 'Group'
	- src sample directory: contanining a Homer Tag directory, TSV, i.e., <src sample dir>/<sampleName>/TSV
	- des sample directory: destination sample directory to put sample (i.e. group) folders, i.e., <des sample dir>/<groupName>/TSV
Options:
	-t: if set, dry run simply displaying pooling message, default=off
	-b: if set, bsub are submitted for jobs, default=off
	-r : If set, perform robust estimation of fragment length, default=OFF
		because Homer gives unreliable fragment length when fragment length is close to the read length" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
testOnly=FALSE
bsub=FALSE
robust=FALSE
while getopts ":brt" opt; do
	case $opt in
		t)
			testOnly=TRUE
			;;
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


echo -e "Pooling replicate bam files" >&2
echo -e "  - sampleInfo: $sampleInfo" >&2
echo -e "  - srcDir: $srcDir" >&2
echo -e "  - desDir: $desDir" >&2
echo -e "  - testOnly: $testOnly" >&2
echo -e "  - bsub:   $bsub" >&2
echo -e "  - robust: $robust" >&2
echo -e "" >&2

groupL=`tail -n +2 $sampleInfo | grep -v -e ^$ -e ^# -e ^Id | cut -f 3 | sort | uniq`

mkdir -p $desDir
log=${desDir}/pool.log
echo -ne "" > $log

for group in ${groupL[@]}
do
	echo -e "Processin $group" >&2
	des=${desDir}/${group}/TSV

	## List of replicate tag directories
	srcL=( `tail -n +2 $sampleInfo | grep -v -e ^$ -e "^#" | gawk '{ if($3 == "'$group'") printf "'$srcDir'/%s/TSV\n", $2 }'` )
	assertDirExist ${srcL[@]}

	if [ -d $des ] && [ -f ${des}/Autocorrelation.png ];then
		echo -e "  - Warning: $des already exists, skip" >&2
		continue
	fi


	## Merging to destination
	echo -e "  - Creating $des" >&2

	if [ "$testOnly" = "TRUE" ];then
		for src in ${srcL[@]}
		do
			echo -e "\t$src" >&2
		done 
		continue
	else
		mkdir -p $desDir
		echo -e "$group" >> $log
		for src in ${srcL[@]}
		do
			echo -e "\t$src" >&2
		done 2>&1 | tee -a $log
	fi

	## info.txt creation
	sampleDir=`dirname $des`
	mkdir -p $sampleDir
	echo -e "  - Creating ${sampleDir}/info.txt" >&2
	echo ${group} > ${sampleDir}/info.txt

	## Pooling
	if [ ${#srcL[@]} -eq 1 ];then
		echo -e "  - Only single sample; Creating a symbolic link" >&2
		ln -vsr ${srcL[@]} $des;
	else
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
	fi
done

