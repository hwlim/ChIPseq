#!/usr/bin/env bash


source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT


function printUsage {
	echo "Usage: `basename $0` (options) <tagDir>
Description: run peak calling using homer in factor mode using a given tag directory
Output:
	- <outPrefix>.homer.txt
	- <outPrefix>.homer.bed
	- <outPrefix>.homer.exBL.bed
	- <outPrefix>.homer.exBL.1rpm.bed
	- <outPrefix>.homer.log
Options:
	-o <outPrefix>: output file prefix including path. default=homerPeak
	-m <mask>: blacklist bed file (optional)
	-i <input tag dir>: input (control) tag directory, default=NULL
	-s <optStr>: comma-separated findPeak option string. ex) -fragLength,150,-size,200. default=NULL" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
outPrefix=""
blacklist="NULL"
input="NULL"
optStr="NULL"
while getopts ":o:m:i:s:" opt; do
	case $opt in
		o)
			outPrefix=$OPTARG
			;;
		m)
			blacklist=$OPTARG
			;;
		i)
			input=$OPTARG
			;;
		s)
			optStr=$OPTARG
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

chip=$1
assertDirExist $chip

if [ "$blacklist" != "NULL" ];then
	assertFileExist $blacklist
fi
if [ "$input" != "NULL" ]; then
	assertDirExist $input
fi

optStr=`echo $optStr | sed 's/,/ /g'`

#############################################
# Output files
desRaw=${outPrefix}.homer.txt
desBed=${outPrefix}.homer.bed
desBedexBL=${outPrefix}.homer.exBL.bed
desBedexBL1rpm=${outPrefix}.homer.exBL.1rpm.bed
log=${outPrefix}.homer.log

#############################################
# Peak calling
echo -e "##########################################" >&2
echo -e "Running Homer findPeaks: Factor mode" >&2
echo -e "  outPrefix = $outPrefix" >&2
echo -e "  ChIP = $chip" >&2
echo -e "  Input = $input" >&2
echo -e "  Blacklist = $blacklist" >&2
echo -e "  optStr = \"$optStr\"" >&2
echo -e "" >&2

if [ "$optStr" = NULL ];then
	optStr=""
fi
desDir=`dirname $outPrefix`
mkdir -p $desDir

if [ "$input" = NULL ];then
	echo -e "  1) findPeaks $chip -o $desRaw $optStr -style factor" >&2
	findPeaks $chip -o ${desRaw} $optStr -style factor 2>&1 | tee ${log}
else
	echo -e "  1) findPeaks $chip -o ${desRaw} -i $input $optStr -style factor" >&2
	findPeaks $chip -o ${desRaw} -i $input $optStr -style factor 2>&1 | tee ${log}
fi

echo -e "\tConverting to a bed file" >&2
cat $desRaw \
	| grep -v -e ^# -e ^$ \
	| gawk '{ printf "%s\t%d\t%d\t%s\t%s\t%s\n", $2,$3,$4,$1,$6,$5 }' \
	> $desBed

#############################################
# Peak post-processing
if [ "$blacklist" != "NULL" ] && [ -f "$blacklist" ];then
	echo -e "  2) Filtering ENCODE blacklist regions" >&2
	echo -e  "intersectBed -a $desBed -b $blacklist -v > $desBedexBL" >&2
	intersectBed -a $desBed -b $blacklist -v > $desBedexBL
else
	echo -e "No black list file exists, pass" >&2
	cp $desBed $desBedexBL
fi

cat $desBedexBL | gawk '$5 > 1' > $desBedexBL1rpm
N=`cat $desBedexBL1rpm | wc -l`
echo -e "  >> Final number of peaks (>1RPM/exBL): $N" >&2
echo -e "" >&2
