#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT


function printUsage {
	echo "Usage: `basename $0` (options) " >&2
	echo "Description: run peak calling using homer in factor mode using a given tag directory" >&2
	echo -e "Options:" >&2
	echo -e "\t-o <desDir>: destination directory (mandatory)" >&2
	echo -e "\t-g <genome>: genome (mandatory)" >&2
	echo -e "\t-i <input tag dir>: input (control) tag directory, default=NULL" >&2
	echo -e "\t-n <name>: name, <desDir>/<name>.homer becomes a common prefix for output files. default=homerPeak" >&2
	echo -e "\t-s <optStr>: comma-separated findPeak option string. ex) -fragLength,150,-size,200. default=NULL" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
desDir=""
genome=""
ctrl="NULL"
name="homerPeak"
optStr=""
while getopts ":o:g:s:i:n:" opt; do
	case $opt in
		o)
			desDir=$OPTARG
			;;
		g)
			genome=$OPTARG
			;;
		s)
			optStr=$OPTARG
			;;
		i)
			ctrl=$OPTARG
			;;
		n)
			name=$OPTARG
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


#######################################################
## Option verification
if [ "$genome" = "" ];then
	echo "Error: genome must be specified (-g)" >&2
	printUsage
	exit 1
fi

if [ "$ctrl" != "NULL" ];then
	assertDirExist $ctrl
fi

if [ "$optStr" = NULL ];then
	optStr="-norm 1000000 "
else
	optStr=`echo $optStr | sed 's/,/ /g'`
	optStr="$optStr -norm 1000000"
fi

chip=$1
assertDirExist $chip

blacklist=~/Research/Common_Data/${genome}/ENCODE-blacklist.bed
assertFileExist $blacklist


#############################################
## Main code for peak calling


## Output files
des=${desDir}/${name}.homer.txt
desBed=${desDir}/${name}.homer.bed
desBed_exBL=${desDir}/${name}.homer.exBL.bed
desBed_exBL_1rpm=${desDir}/${name}.homer.exBL.1rpm.bed
log=${desDir}/${name}.homer.log

## Start
mkdir -p $desDir

echo -e "##########################################" >&2
echo -e "Running Homer findPeaks: Factor mode" >&2
echo -e "  Name = $name" >&2
echo -e "  ChIP = $chip" >&2
echo -e "  Input = $ctrl" >&2
echo -e "  desDir = $desDir" >&2
echo -e "  optStr = $optStr" >&2
echo -e "" >&2

echo "##########################################
Running Homer findPeaks: Factor mode
  Name = $name
  ChIP = $chip
  Input = $ctrl
  desDir = $desDir
  optStr = $optStr" > $log
echo "" >> $log
echo "" >> $log


tmp=__temp__.$$.txt
## Peak calling
if [ "$ctrl" = "NULL" ];then
	echo -e "  1) findPeaks $chip -o ${desDir}/${name}.homer.txt $optStr -style factor" >&2
	findPeaks $chip -o ${tmp} $optStr -style factor 2>&1 | tee -a ${log}
else
	echo -e "  1) findPeaks $chip -o ${des} -i $ctrl $optStr -style factor" >&2
	findPeaks $chip -o ${tmp} -i $ctrl $optStr -style factor 2>&1 | tee -a ${log}
fi
mv $tmp $des
#chrM-1	chrM	15946	16146	+	27262.1	0.778	156573.000000	15.50	0.00e+00	0.61

cat $des \
	| grep -v -e ^# -e ^$ \
	| gawk '{ printf "%s\t%d\t%d\t%s\t%s\t+\n", $2,$3,$4,$1,$6 }' \
	> $desBed
#pos2bed.pl $des | grep -v ^# | grep -v ^$ > $desBed


## Peak post-processing
#echo -e "  2) Normalized tag counting (RPM)" >&2
#countTagsBed2.sh -q -ra $desBed $chip \
#	| gawk '{ printf "%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,$8,$6 }' \
#	| sort -k5,5nr \
#	>
#	> $desBed

echo -e "  2) Filtering ENCODE blacklist regions" >&2
if [ -f $blacklist ];then
	intersectBed -a $desBed -b $blacklist -v > $desBed_exBL
else
	echo -e "No black list file exists, pass" >&2
	cp $desBed $desBed_exBL
fi

cat $desBed_exBL | gawk '$5 > 1' > $desBed_exBL_1rpm
N=`cat $desBed_exBL_1rpm | wc -l`
echo -e "  >> Final number of peaks (>1RPM/exBL): $N" >&2

