#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT

function printUsage {
	echo "Usage: `basename $0` (options) [alignFile]
Description: Make data and homer tag directory using a given alignment file
Options:
	-o : dataDir, default=<src file name>
	-n : name, default=<src file name>
	-t : -tbp option for homer, default=0
	-c : regular expression of chr to select, default=NULL
		e.g. 'chr[0-9XY]*$' to select autosomes and sex chromosomes
Input:
	- bed file, plain or gzipped
	- bam
	- sam
Output:
	<dataDir>/TSV<tbp>: tagDir
	<dataDir>/info.txt: information file containing <name>">&2
}

###################################
## option and input file handling
dataDir=""
name=""
tbp=0
chrRegex=NULL
while getopts ":o:n:t:c:" opt; do
	case $opt in
		o)
			dataDir=$OPTARG
			;;
		n)
			name=$OPTARG
			;;
		t)
			tbp=$OPTARG
			;;
		c)
			chrRegex=$OPTARG
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
assertFileExist $src

tagDir=${dataDir}/TSV${tbp}
srcFile=`basename $src`

if [ "$dataDir" = "" ];then
	echo -e "Error: dataDir (-o) must be specified" >&2
	exit 1
fi

if [ "$name" = "" ];then
	name="$srcFile"
fi

if [ ! ${tbp} -ge 0 ];then
	echo -e "Error: -t must be >= 0" >&2
	exit 1
fi

echo -e "============================================"
echo -e "Making Homer Tag Directories"
echo -e "  src: $src" >&2
echo -e "  dataDir: $dataDir" >&2
echo -e "  tagDir: $tagDir" >&2
echo -e "  tbp: $tbp" >&2
echo -e "  chrRegex: $chrRegex" >&2
echo -e "" >&2

printAlign(){
	local src=$1
	local ext=${src##*.}

	if [ "$chrRegex" = "NULL" ];then
		if [ "${ext}" = "bed" ];then
			cat $src
		elif [ "${ext}" = "gz" ];then
			zcat $src
		elif [ "${ext}" = "bam" ];then
			bamToBed -i $src
		elif [ "${ext}" = "sam" ];then
			cat $src
		else
			echo "Error: Unknown file format: ${1}" >&2
			exit 1
		fi
	else
		if [ "${ext}" = "bed" ];then
			cat $src | gawk '$1 ~ /'$chrRegex'/'
		elif [ "${ext}" = "gz" ];then
			zcat $src | gawk '$1 ~ /'$chrRegex'/'
		elif [ "${ext}" = "bam" ];then
			bamToBed -i $src | gawk '$1 ~ /'$chrRegex'/'
		elif [ "${ext}" = "sam" ];then
			cat $src | gawk '$3 ~ /'$chrRegex'/'
		else
			echo "Error: Unknown file format: ${1}" >&2
			exit 1
		fi
	fi
}


ext=${src##*.}
if [ "${ext}" = "sam" ];then
	optStr="-format sam"
elif [ "${ext}" = "bed" ] || [ "${ext}" = "gz" ] || [ "${ext}" = "bam" ]; then
	optStr="-format bed"
else
	echo -e "Error: invalid file format with an extension: $ext" >&2
	exit 1
fi

if [ $tbp -gt 0 ];then
	optStr="$optStr -tbp $tbp"
fi

echo "1) Creating tag directory: ${tagDir}" >&2
mkdir -p $tagDir
echo "$name" > ${dataDir}/info.txt

echo -e "printAlign $src | makeTagDirectory ${tagDir} /dev/stdin ${optStr}" >&2
( printAlign $src \
	| makeTagDirectory ${tagDir} /dev/stdin ${optStr} ) 2>&1 \
	| tee ${tagDir}/TSV.log

drawAutoCorrplot.r -t ${name} ${tagDir}


