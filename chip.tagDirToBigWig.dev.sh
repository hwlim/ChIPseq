#!/usr/bin/env bash

trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT
source $COMMON_LIB_BASE/commonBash.sh

function printUsage {
	#echo "***WATNING: chrM coordination handling must be further improved***"
	echo -e "Usage: `basename $0` (options) [tagDir]
Options:
	-o <out>: Output bigWig file with an extension bw or bigwig, required
	-g <chrom>: chromosome size file. must be sorted in compliance with bedtool" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


des=""
chrom=""
sortMem=5G
while getopts ":o:g:m:" opt; do
	case $opt in
		o)
			des=$OPTARG
			;;
		g)
			chrom=$OPTARG
			;;
		m)
			sortMem=$OPTARG
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
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

if [ "$des" = "" ];then
	echo "Error: destimation file must be specified" >&2
	printUsage
	exit 1
fi
if [ "$chrom" = "" ];then
	echo "Error: chrom.size file must be specified" >&2
	printUsage
	exit 1
fi

tagDir=$1
assertFileExist $tagDir

tmpBg=${des%.bw}.bg

echo -e "===============================================" >&2
echo -e "Creating bigWig file from Homer tag directory" >&2
echo -e "  - src = $tagDir" >&2
echo -e "  - des = $des" >&2
echo -e "  - chrom = $chrom" >&2

echo -e "1) Making bedGraph file" >&2
makeUCSCfile $tagDir -strand both -norm 1000000 -fsize 1e10 \
	| grep "^chr" \
	| gawk 'BEGIN{
			split("'$chrLen'", lenAr, ",");
			for(i=1;i<length(lenAr);i=i+2){
				lenD[lenAr[i]]=lenAr[i+1];
			}
		}
		{
			if($3>lenD[$1] || $2<0){
				next
			}else{
				print $0
			}
		}'\
	| gawk '{printf "%s\t%s\t%s\t%.5f\n", $1,$2,$3,$4}'\
	| sort -S $sortMem -k1,1 -k2,2n -k3,3nr \
	> ${tmpBg}

echo -e "2) Converting to bigWig file" >&2
bedGraphToBigWig ${tmpBg} ${chrom} ${des}
rm ${tmpBg}


