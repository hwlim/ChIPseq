#!/usr/bin/env bash

# Sort chromosome order of given bam file(s) lexicographically

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT

function printUsage {
	echo "Usage: `basename $0` (options) [srcBam] [desBam]" >&2
	echo -e "Options:" >&2
        echo -e "\t-m <memory>: memory size for sorting in K/M/G. default=2G" >&2
#        echo -e "\t-f : flag option" >&2
}



###################################
## option and input file handling
mem=2G
while getopts ":m:" opt; do
	case $opt in
		m)
			mem=$OPTARG
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
if [ $# -lt 2 ];then
	printUsage
	exit 1
fi

src=$1
des=$2

echo -e "==============================================" >&2
echo -e "Sorting bam file chromosomes lexicographically" >&2
echo -e "  src = $src" >&2
echo -e "  des = $des" >&2
echo -e "  mem = $mem" >&2

echo -e "  1) Sorting " >&2
printBam(){
	samtools view -H $1 | grep "^@HD"
	samtools view -H $1 | grep "^@SQ" | sort -k2,2
	samtools view -H $1 | grep -v -e "^@SQ" -e "^@HD"
	samtools view $1
}

printBam $src \
	| samtools view -b - \
	| samtools sort -m $mem -O BAM -T __temp__.$$ - \
	> __temp__.$$.bam
mv __temp__.$$.bam $des

echo -e "  2) Indexing " >&2
samtools index $des
