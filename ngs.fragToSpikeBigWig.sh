#!/usr/bin/env bash


source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [bed]
Description:
	Make a spikein BigWig file from a fragment BED file in RPM scale (default) or manually scaled
Input:
	Fragment file that contains both target and spikein fragments and
	spikein fragments are prefixed with a specific keyword such as dm-
Output:
	RPM-normalizzed bigwig file with prefix deleted in the chromosome
Options:
	-o <outFile>: Destination file. required
	-l <minLen>: lower bound of fragment length to use (including the boundary), default=0
	-L <maxLen>: upper bound of fragment length to use (including the boundary), default=1000000 (infinite)
	-r <resize>: resize to the give length around the center, default=-1 (no resize)
	-g <chromSize>: spikein chromosome size file without prefix
	-p <spikeinPrefix>: spikein prefix. required
	-m <memory>: memory size for sorting bedGraph file, default=5G
	-s <scale factor>: Manual scaling factor. This value is multiplied to \"raw read count\" primarily for spike-in based scaling.
			If 0, RPM normalized. default=0" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
des=NULL
minLen=0
maxLen=1000000
resize=-1
genome=NULL
memory=5G
scaleFactor=0
spikeinPrefix=.
while getopts ":o:l:L:r:g:m:s:p:" opt; do
	case $opt in
		o)
			des=$OPTARG
			;;
		l)
			minLen=$OPTARG
			;;
		L)
			maxLen=$OPTARG
			;;
		r)
			resize=$OPTARG
			;;
		g)
			genome=$OPTARG
			;;
		m)
			memory=$OPTARG
			;;
		s)
			scaleFactor=$OPTARG
			;;
		p)
			spikeinPrefix=$OPTARG
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

if [ "$genome" = "NULL" ];then
	echo -e "Error: chromosome size (-g) must be specified" >&2
	exit 1
fi
src=$1
assertFileExist $src

if [ "$des" == "NULL" ];then
	echo -e "Error: Destination file (-o) must be specified" >&2
	exit 1
fi

assertFileExist $genome

###################################
## main code
printBed(){
	local src=$1
	assertFileExist $src

	local ext=${src##*.}
	
	# avoid using zero coordinate reads
	if [ "$ext" == "gz" ];then
		zcat $src | gawk '$2 > 0'
	else
		cat $src | gawk '$2 > 0'
	fi
}

printFrag(){
	local src=$1
	assertFileExist $src

	printBed $src \
		| gawk 'BEGIN{
				minLen='$minLen';
				maxLen='$maxLen';
				resize='$resize';
			} $1 ~ /^'$spikeinPrefix'/ {
				fragLen=$3-$2
				if( fragLen < minLen || fragLen > maxLen ) next

				if( resize > -1 ){
					c = ($3+$2)/2
					h = resize / 2
					s = c - h
					e = c + h
					if(s < 0) s = 0
					printf "%s\t%d\t%d\tNULL\t0\t%s\n", $1, s, e, $6
				}else{
					printf "%s\t%d\t%d\tNULL\t0\t%s\n", $1, $2, $3, $6
				}
			}' \
		| sed 's/^'$spikeinPrefix'//'
}

desDir=`dirname $des`
mkdir -p $desDir

tmpBG=${TMPDIR}/__temp__.$$.bedGraph
tmpBW=${TMPDIR}/__temp__.$$.bw

echo -e "Creating BigWig file from a fragment bed file
  - src = $src
  - des = $des
  - fragLen = $minLen - $maxLen (bp)
  - resize = $resize (bp)
  - chromSize = $genome
  - spikeinPrefix = $spikeinPrefix" >&2

if [ $scaleFactor == "0" ];then
	echo -e "  1) Calculating scale factor for RPM normalization" >&2
	ttc=`printFrag $src | wc -l`
	scaleFactor=`echo $ttc | gawk '{ printf "%f", 1000000/$1}'`
	echo -e "\tTTC = $ttc (scaleFactor $scaleFactor)" >&2
else
	echo -e "  1) Scale factor was manually assigned" >&2
	echo -e "\tscaleFactor: ${scaleFactor}" >&2
fi

echo -e "  2) Making bedGraph file" >&2
printFrag $src \
	| sort -S $memory -k1,1 -k2,2n -k3,3n \
	| genomeCoverageBed -bg -scale $scaleFactor -g $genome -i stdin \
	| gawk '{ printf "%s\t%s\t%s\t%.5f\n", $1,$2,$3,$4 }' \
	> $tmpBG

echo -e "  3) Converting to bigWig file" >&2
bedGraphToBigWig ${tmpBG} $genome ${tmpBW}
mv ${tmpBW} ${des}

