#!/usr/bin/env bash

###########################################3
# Cut&Run tools
# Written by Hee-Wooong Lim
# 
# Pool bam files of multiple replicates by group
#
# To do:
# - Handling premature stop by deleting incomplete output file

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) <sample.tsv> <src sample directory> <des sample directory>
Description:
	Merge multiple fragment bed files of a group according to sample/group information within a given sample.tsv file
Input:
	- sample.tsv file: containing columns 'Name' and 'Group'
	- src fragment directory: contanin replicate fragment bed files
	  e.g. <sample dir>
	  ├── <sample1>.frag.bed.gz
	  └── <sample2>.frag.bed.gz
	  In default, input fragments files are assumed to be sorted. So merge-sort will be performed to create output.
	- des fragment directory: to write merged fragment bed files
	  e.g. <group dir>
	  ├── <group1>.frag.bed.gz
	  └── <group2>.frag.bed.gz

Options:
	-t: if set, dry run simply displaying pooling message, default=off
	-b: if set, bsub are submitted for merging bam files, default=off
	-u: if set, input/output fragment files are assumed to be 'unsorted', not merge sort performed. default=off
	-f: if set, force overwrite existing destination bam files, default=off" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
bsub=FALSE
testRun=FALSE
unsorted=FALSE
overwrite=FALSE
while getopts ":buft" opt; do
	case $opt in
		t)
			testRun=TRUE
			;;
		b)
			bsub=TRUE
			;;
		u)
			unsorted=TRUE
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

groupL=`cat $sampleInfo | grep -v -e ^$ -e "^#" -e ^Id | cut -f 3 | sort | uniq`

echo -e "Pooling replicate fragment bed files" >&2
echo -e "  - sampleInfo: $sampleInfo" >&2
echo -e "  - srcDir: $srcDir" >&2
echo -e "  - desDir: $desDir" >&2
echo -e "  - unsorted: $unsorted" >&2
echo -e "  - bsub:   $bsub" >&2

mkdir -p $desDir
for group in ${groupL[@]}
do
	des=${desDir}/${group}/fragment.bed.gz
	log=${desDir}/${group}/fragment.log

	## Checking existing destination file
	if [ -f $des ] && [ "$overwrite" != "TRUE" ];then
		echo -e "Warning: $des already exists, pass" >&2
		continue
	fi

	## List of replicate bam files
	srcL=( `cat $sampleInfo | grep -v -e ^$ -e "^#" -e ^Id | gawk '{ if($3 == "'$group'") printf "'$srcDir'/%s/fragment.bed.gz\n", $2 }'` )
	assertFileExist ${srcL[@]}

	## Merging to destination
	echo -e "Creating $des" >&2
	for src in ${srcL[@]}
	do
		echo -e "  - $src" >&2
	done

	if [ "$testRun" == "TRUE" ];then
		continue
	fi

	mkdir -p ${desDir}/${group}

	echo -ne "" > $log
	for src in ${srcL[@]}
	do
		echo -e "- $src" >> $log
	done

	# ## command to uncompress
	# cmd="sort -m -k1,1 -k2,2n -k3,3n"
	# for src in ${srcL[@]}
	# do
	# 	cmd="${cmd} <( zcat $src )"
	# done
	# echo -e "Submitting" >&2
	# echo "  >> $cmd" >&2

	# if [ "$bsub" == "TRUE" ];then
	# 	## Parallel processing using HPC:lsf
	# 	bsub -W 24:00 -n 1 "eval $cmd | gzip > ${TMPDIR}/__temp__.$$; mv ${TMPDIR}/__temp__.$$ $des"
	# else
	# 	## Sequential processing
	# 	eval $cmd | gzip > ${TMPDIR}/__temp__.$$
	# 	mv ${TMPDIR}/__temp__.$$ $des
	# 	#eval $cmd | gzip > $des
	# 	#ngs.concateBamFiles.sh -o $des ${srcL[@]}
	# fi
	# echo "" >&2

	tmp=${TMPDIR}/__temp__.$$.${group}.bed.gz
	if [ "$unsorted" == "TRUE" ];then
		if [ "$bsub" == "TRUE" ];then
			## Parallel processing using HPC:lsf
			bsub -W 24:00 -n 1 "cat ${srcL[@]} > $tmp; mv $tmp $des"
		else
			## Sequential processing
			cat ${srcL[@]} > $tmp
			mv $tmp $des
		fi
	else
		inputStr=""
		for src in ${srcL[@]}
		do
			inputStr="${inputStr} <( zcat $src )"
		done

		if [ "$bsub" == "TRUE" ];then
			## Parallel processing using HPC:lsf
			bsub -W 24:00 -n 1  <<- EOF
				#!/usr/bin/env bash
				sort -m -k1,1 -k2,2n -k3,3n ${inputStr} | gzip > $tmp
				mv $tmp $des
			EOF
		else
			## Sequential processing
			#echo -e "sort -m -k1,1 -k2,2n -k3,3n ${inputStr}"
			eval "sort -m -k1,1 -k2,2n -k3,3n ${inputStr} | gzip > $tmp"
			mv $tmp $des
		fi
	fi
	echo "" >&2

done

