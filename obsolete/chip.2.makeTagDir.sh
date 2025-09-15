#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh

source ./config.ChIP.star.sh

## data list file
if [ $# -gt 0 ];then
	srcDataList=$1
	assertFileExist $1
else
	echo -e "Error: Data list file must be specified" >&2
	exit 1
fi

if [ "${alignDir2}" != NULL ];then
	srcBase=${alignDir2}
else
	srcBase=${alignDir}
fi
desBase=.
isDirExist $srcBase
mkdir -p $desBase

echo -e "============================================"
echo -e "Making Homer Tag Directories"
echo -e "  SrcBase: $srcBase" >&2
echo -e "  DesBase: $desBase" >&2
echo -e "  Genome: $genome" >&2
echo -e "  chrRegex: $chrRegex" >&2
echo -e "============================================"
echo -e "" >&2

printAlign(){
	local src=$1
	local ext=${src##*.}

	if [ -z ${chrRegex+x} ];then
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


while IFS='' read -r line || [[ -n "$line" ]];
do
	[[ "$line" =~ ^#.*$ ]] &&  continue
	[[ "$line" = "" ]] && continue
	N_field=`echo "$line" | gawk -F "\t" '{ print NF }'`
	if [ "$N_field" -eq 5 ];then
		src=${srcBase}/${id}.bam
	elif [ "$N_field" -eq 6 ];then
		alignFile=`echo "$line" | cut -d "	" -f 6`
		if [ "${alignFile}" = "NULL" ];then
			src=${srcBase}/${id}.bam
		else
			src=${srcBase}/${alignFile}
		fi
	else
		echo -e "Error: Invalid number of columns, $N_field" >&2
		exit 1
	fi
	
	id=`echo "$line" | cut -d "	" -f 1`
	dataDir=`echo "$line" | cut -d "	" -f 2`
	name=`echo "$line" | cut -d "	" -f 3`

	
	
	echo -e "Processing $src" >&2
	echo -e "  id: ${id}" >&2
	echo -e "  name: ${name}" >&2
	echo -e "  dataDir: ${dataDir}" >&2
	echo -e "  src: ${src}" >&2
	echo -e "" >&2

	## name 
	echo "0) Creating ${dataDir}/info.txt" >&2
	mkdir -p ${dataDir}
	echo "$name" > ${dataDir}/info.txt


	ext=${src##*.}
	if [ "${ext}" = "sam" ];then
		optStr="-format sam"
	else
		optStr="-format bed"
	fi
	
	echo "1) Creating tag directory: ${dataDir}/TSVrd" >&2
	mkdir -p ${dataDir}/TSVrd

	echo -e "\t printAlign $src | makeTagDirectory ${dataDir}/TSVrd /dev/stdin ${optStr}" >&2
	( printAlign $src \
		| makeTagDirectory ${dataDir}/TSVrd /dev/stdin ${optStr} ) 2>&1 \
		| tee ${dataDir}/TSVrd/TSVrd.log
	drawHomerAutoCorr.r -t ${name} ${dataDir}/TSVrd

	echo "2) Creating tag directory: ${dataDir}/TSVu" >&2
	mkdir -p ${dataDir}/TSVu
	echo -e "tagDir2bed.pl ${dataDir}/TSVrd | makeTagDirectory ${dataDir}/TSVu /dev/stdin -format bed" >&2
	tagDir2bed.pl ${dataDir}/TSVrd \
		| makeTagDirectory ${dataDir}/TSVu /dev/stdin -format bed 2>&1 | tee ${dataDir}/TSVu/TSVu.log
	drawHomerAutoCorr.r -t ${name} ${dataDir}/TSVu

	echo -e "" >&2
	
done < $srcDataList

