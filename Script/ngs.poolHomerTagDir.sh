#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT

function printUsage {
	echo "Usage: `basename $0` (options) [tagDir1] [tagDir2] ...
Description: Create a homer tag directory using a given alignment file
Input:
	A list of Homer tag directories
Options:
	-o : output tag directory,required
	-t : -tbp option for homer, default=0
	-r : If set, perform robust estimation of fragment length, default=OFF
	     because Homer gives unreliable fragment length when the max-autocorr length is too small (e.g. CUT&RUN)" >&2
}

###################################
## option and input file handling
outiDir=NULL
tbp=0
robust=FALSE
while getopts ":o:t:r" opt; do
	case $opt in
		o)
			outDir=$OPTARG
			;;
		t)
			tbp=$OPTARG
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


srcL=( $@ )
assertDirExist ${srcL[@]}

#outDir=${dataDir}/TSV${tbp}
#srcFile=`basename $src`

if [ "$outDir" = "NULL" ];then
	echo -e "Error: outDir (-o) must be specified" >&2
	exit 1
fi

if [ ! ${tbp} -ge 0 ];then
	echo -e "Error: -t must be >= 0" >&2
	exit 1
fi

echo -e "============================================"
echo -e "Pooling Homer Tag Directories"
echo -e "  src:" >&2
for src in ${srcL[@]}
do
	echo -e "    - $src" >&2
done
echo -e "  outDir: $outDir" >&2
echo -e "  tbp: $tbp" >&2
echo -e "  robust: $robust" >&2
echo -e "" >&2


if [ $tbp -gt 0 ];then
	optStr="-tbp $tbp"
else
	optStr=""
fi

mkdir -p $outDir

( makeTagDirectory ${outDir} -d ${srcL[@]} ${optStr} ) 2>&1 \
	| tee ${outDir}/TSV.log

# drawHomerAutoCorr.r -t ${name} ${outDir} # not used because "name" is not known here


if [ "$robust" == "TRUE" ];then
	echo -e "Performing robust fragment estimation using Homer tagAutocorrelation.txt" >&2
	srcAuto=${outDir}/tagAutocorrelation.txt
	assertFileExist $srcAuto
	fragLen=`estimateHomerFragLen.r $srcAuto`
	makeTagDirectory ${outDir} -update -fragLength $fragLen 2>&1 | tee ${outDir}/TSV.log
fi
