#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT

function printUsage {
	echo "Usage: `basename $0` (options) [alignFile]
Description: Create a homer tag directory using a given alignment file
Input:
	BAM/SAM/BED/BED.gz
Options:
	-o : output tag directory,required
	-t : -tbp option for homer, default=0
	-c : regular expression of chr to select, default=NULL
		e.g. 'chr[0-9XY]*$' to select autosomes and sex chromosomes
	-r : If set, perform robust estimation of fragment length, default=OFF
	     because Homer gives unreliable fragment length when the max-autocorr length is too small (e.g. CUT&RUN)" >&2
}

###################################
## option and input file handling
outDir=NULL
tbp=0
chrRegex=NULL
robust=FALSE
while getopts ":o:t:c:r" opt; do
	case $opt in
		o)
			outDir=$OPTARG
			;;
		t)
			tbp=$OPTARG
			;;
		c)
			chrRegex=$OPTARG
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
if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


src=$1
assertFileExist $src

#outDir=${dataDir}/TSV${tbp}
srcFile=`basename $src`

if [ "$outDir" = "NULL" ];then
	echo -e "Error: outDir (-o) must be specified" >&2
	exit 1
fi

if [ ! ${tbp} -ge 0 ];then
	echo -e "Error: -t must be >= 0" >&2
	exit 1
fi

echo -e "============================================"
echo -e "Making Homer Tag Directories"
echo -e "  src: $src" >&2
echo -e "  outDir: $outDir" >&2
echo -e "  tbp: $tbp" >&2
echo -e "  chrRegex: $chrRegex" >&2
echo -e "  robust: $robust" >&2
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

echo "1) Creating tag directory: ${outDir}" >&2
mkdir -p $outDir

echo -e "printAlign $src | makeTagDirectory ${outDir} /dev/stdin ${optStr}" >&2
( printAlign $src \
	| makeTagDirectory ${outDir} /dev/stdin ${optStr} ) 2>&1 \
	| tee ${outDir}/TSV.log

#drawHomerAutoCorr.r -t ${name} ${outDir}


if [ "$robust" == "TRUE" ];then
	echo -e "Performing robust fragment estimation using Homer tagAutocorrelation.txt" >&2
	srcAuto=${outDir}/tagAutocorrelation.txt
	assertFileExist $srcAuto
	fragLen=`estimateHomerFragLen.r $srcAuto`
	makeTagDirectory ${outDir} -update -fragLength $fragLen 2>&1 | tee ${outDir}/TSV.log
fi
