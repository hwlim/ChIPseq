#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh

function printUsage {
	echo -e "Process single alginFile: tag directory, peak calling, ucsc bedGraph" >&2
	echo -e "Input format: [chr] [start] [end] [name] [count] [direc]" >&2
	echo -e "\tSAM/BAM/BED/bed.gz file" >&2
	echo -e "Usage: dataPipeline2.sh (options) [alignFile]" >&2
	echo -e "Options:" >&2
	echo -e "\t-n <name>: data name" >&2
	echo -e "\t-g <genome>: genome" >&2
	echo -e "\t-o <desDir>: destination directory" >&2
	echo -e "\t-1: alignment file is 1-based, not 0-based like BED format.\n\t\tIf set, convert it into zero-based before processing" >&2
	echo -e "\t-t: tag directory creation, TSVrd & TSVu" >&2
	echo -e "\t-p: basic peak calling without input" >&2
	echo -e "\t-b: generate bed graph" >&2
#	echo -e "\t-w: generate bigWig" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi

name=""
genome=""
desDir=.
toTSV=0
toPeak=0
toBedGraph=0
toBigWig=0
oneBase=0
while getopts ":n:g:o:1tpbw" opt; do
	case $opt in
		n)
			name=$OPTARG
			;;
		g)
			genome=$OPTARG
			;;
		o)
			desDir=$OPTARG
			;;
		1)
			oneBase=1
			;;
		t)
			toTSV=1
			;;
		p)
			toPeak=1
			;;
		b)
			toBedGraph=1
			;;
		w)
			toBigWig=1
			;;
		/?)
			echo "Invalid options: -$OPTARG" >&2
			printUsage
			exit 1
			;;
		:)
			echo "Options -$OPTARG requires an argument." >&2
			printUsage
			exit 1
			;;
	esac
done

#echo toPeak=$toPeak
shift $((OPTIND-1))

if [ "$genome" == "" ];then
	echo "Error: genome is required." >&2
	printUsage
	exit 1
fi

if [ $# -lt 1 ];then
	echo "Error: at least one alignFile is required.">&2
	printUsage
	exit 1
fi

if [ $toTSV -eq 1 ];then
	assertFileExist $1
fi
alignFile=$1

# Option print
echo -e "Align file:\t${alignFile}" >&2
echo -e "Name:\t\t${name}" >&2
echo -e "Genome:\t\t${genome}">&2
echo -e "Destination:\t${desDir}">&2
echo -e "Make TagDir:\t${toTSV}">&2
echo -e "Do Peak call:\t${toPeak}">&2
echo -e "Make BedGraph:\t${toBedGraph}">&2

mkdir -p ${desDir}

if [ "$name" != "" ];then
	echo "Creating ${desDir}/info.txt" >&2
	echo "$name" > ${desDir}/info.txt
else
	if [ -f ${desDir}/info.txt ];then
		name=`cat ${desDir}/info.txt`
	else
		echo "Error: No information on dataName, -n <name> options is required." >&2
		exit 1
	fi
fi


catCmd(){
	local src=$1
	local ext=${src##*.}

	if [ "${ext}" = "bed" ] || [ "${ext}" = "sam" ];then
		catCmd="cat"
	elif [ "${ext}" = "gz" ];then
		catCmd="zcat"
	elif [ "${ext}" = "bam" ];then
		catCmd="bamToBed -i" 
	else
		echo "Error: Unknown file format: ${alignFile}" >&2
		exit 1
	fi

	if [ $oneBase -eq 1 ];then
		echo -e "Input is 1-based format, converting to 0-base" >&2
		$catCmd $src | gawk -F "\t" '{ printf "%s\t%d\t%d\t%s\t%s\t%s\n", $1,$2-1,$3,$4,$5,$6}'
	else
		echo -e "Input is 0-based formate" >&2
		$catCmd $src
	fi

}

if [ $toTSV -eq 1 ];then
	ext=${alignFile##*.}
	if [ "${ext}" = "bed" ] || [ "${ext}" = "sam" ];then
		optStr="-format bed"
	elif [ "${ext}" = "gz" ];then
		optStr="-format bed"
	elif [ "${ext}" = "bam" ];then
		optStr="-format bed"
	else
		echo "Error: Unknown file format: ${alignFile}" >&2
		exit 1
	fi


	echo "Creating tag directory: ${desDir}/TSVrd" >&2
	mkdir -p ${desDir}/TSVrd

	echo -e "catCmd $alignFile | makeTagDirectory ${desDir}/TSVrd /dev/stdin ${optStr}" >&2
	( catCmd $alignFile \
		| makeTagDirectory ${desDir}/TSVrd /dev/stdin ${optStr} ) 2>&1 \
		| tee ${desDir}/TSVrd/TSVrd.log

	#drawGCplot.sh ${desDir}/TSVrd
	#drawAutoCorrplot.sh ${desDir}/TSVrd
	drawAutoCorrplot.r -t ${name} ${desDir}/TSVrd
	drawTCDistrib.sh ${desDir}/TSVrd
	drawTLDistrib.sh ${desDir}/TSVrd

	echo "Creating tag directory: ${desDir}/TSVu" >&2
	mkdir -p ${desDir}/TSVu
	echo "makeTagDirectory ${desDir}/TSVu -d ${desDir}/TSVrd -tbp 1" >&2
	makeTagDirectory ${desDir}/TSVu -d ${desDir}/TSVrd -tbp 1 2>&1 | tee ${desDir}/TSVu/TSVu.log
	#drawGCplot.sh ${desDir}/TSVu
	#drawAutoCorrplot.sh ${desDir}/TSVu
	drawAutoCorrplot.r -t ${name} ${desDir}/TSVu
	drawTCDistrib.sh ${desDir}/TSVu
	drawTLDistrib.sh ${desDir}/TSVu
fi

if [ $toPeak -eq 1 ];then
	echo "Basic peak calling without input: ${desDir}/HomerPeak_u_NI" >&2
	runHomerPeak.sh -n $name -o ${desDir}/HomerPeak_u_NI ${desDir}/TSVu none
fi

if [ $toBedGraph -eq 1 ];then
	echo "Generating ${name}.bedGraph.gz">&2
	makeBedGraph.sh -o ${name} -g ${genome} -h 40 ${desDir}/TSVu
fi

#mv ${alignFile} ${desDir}
