#!/usr/bin/env bash

###########################################3
# Cut&Run tools
# Written by Hee-Wooong Lim
# 
# Bam file split by target vs spikein

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [bam]
Description:
	Split an input bam file into two: target chromosomes & spikein chromosomes
	Technically, chromosome selection done by ~ operator of gawk
Input:
	Paired-end BAM file
Options:
    -o <outPrefix>: prefix for output files, required
    	<outPrefix>.target.bam
    	<outPrefix>.spikein.bam
    -t <target chromosome regex>: Regular expression for target chromosome selection
		default=^chr[0-9XY]+$ i.e. select all regular and sex chromosomes
    -s <spikein chromosome regex>: Regular expression for spikein chromosome selection
		default=^dm i.e. any chromosome start with 'dm' (note missing $)"
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
outPrefix=""
targetRex="^chr[0-9XY]+$'"
spikeinRex="^dm"
while getopts ":o:t:s:" opt; do
	case $opt in
		o)
			outPrefix=$OPTARG
			;;
		t)
			targetRex=$OPTARG
			;;
		s)
			spikeinRex=$OPTARG
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
src=$1


###################################
## main code

#optStr=""
if [ "$outPrefix" == "" ];then
	echo -e "Error: Output prefix (-o) must be specified" >&2
	exit 1
fi
assertFileExist $src

desTarget=${outPrefix}.target.bam
desSpikein=${outPrefix}.spikein.bam

desDir=`dirname $outPrefix`
mkdir -p $desDir


echo -e "Splitting a BAM file
  - Input: $src
  - Target Chromosome: $targetRex
  - Spikein Chromosome: $spikeinRex
  - outPrefix: $outPrefix
" >&2


tmp=${TMPDIR}/__temp__.$$.bam
#chrList=`samtools view -H $src | sed 's/:/\t/' | gawk '{ if($1=="@SQ" && $2=="SN") print $3 }' | grep -E -w ${chrRegex}`
#samtools view -b -o $tmp $src $chrList 

echo -e "Creating $desTarget" >&2
samtools view -h $src \
	| gawk '$3 ~ /'$targetRex'/ || $1 ~/^@/' \
	| samtools view -b -1 -o $tmp -
mv $tmp $desTarget

echo -e "Creating $desSpikein" >&2
samtools view -h $src \
	| gawk '$3 ~ /'$spikeinRex'/ || $1 ~/^@/' \
	| samtools view -b -1 -o $tmp -
mv $tmp $desSpikein
