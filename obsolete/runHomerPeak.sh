#!/usr/bin/env bash

# Written by Hee Woong Lim, 12/15/2011
#
# Run homer with given options and perform post processing
# Usage:
#	run.homer.sh chipDir inputDir name
#	
#	chipDir/inputDir: chip and input tag directory
#	name: data name, all the results are generated with this name

source $COMMON_LIB_BASE/commonBash.sh

function printUsage {
	echo "Usage: runHomerPeak.sh (options) [chipTagDir] [inputTagDir]" >&2
	echo -e "Options:" >&2
	echo -e "\tGenome independent running" >&2
	echo -e "\tIf no inputTagDir, just enter 'none'." >&2
	echo -e "\t-o <dir>: Output directory relative to cwd (default: cwd)" >&2
	echo -e "\t-r <dir>: Output directory relative to the parent dir of the ChIP tag directory" >&2
	echo -e "\t\t-o option has higher priority than -r." >&2
	echo -e "\t-e <option string enclosed by \">: Extra option string for Homer findPeaks, <default \"-style factor\">" >&2
	echo -e "\t\tin case of histone enriched region, this should be \"-style histone\"" >&2
	echo -e "\t-h: for enrichemnt region finding, ex) histone">&2
	echo -e "\t-n <name>: Out file name" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi

outDirByCWD=""
outDirByTDir=""
optStr="" #-style factor"
isFactor=1
name=""
while getopts ":o:e:n:h" opt;
do
	case $opt in
		o)
			outDirByCWD=$OPTARG
			;;
		r)
			outDirByTDir=$OPTARG
			;;
		e)
			optStr=$OPTARG
			;;
		n)
			name=$OPTARG
			;;
		h)
			isFactor=0
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

if [ $isFactor -eq 1 ];then
	optStr="-style factor $optStr"
else
	optStr="-style histone $optStr"
fi

shift $((OPTIND-1))
if [ $# -ne 2 ];then
	echo "Error: ChIP/Input directories must be given." >&2
	printUsage
	exit 1
fi

if [ "$name" = "" ];then
	echo "Error: -n <name> option must be specified." >&2
	exit 1
fi
chip=$1
isDirExist $chip

input=$2
if [ `echo $input | gawk '{print tolower($1)}'` != "none" ];then
	isDirExist $input
else
	input=none
fi


if [ "$outDirByCWD" = "" ];then
	if [ "$outDirByTDir" = "" ];then
		outDir=.
	else
		outDir=$chip/$outDirByTDir
	fi
else
	outDir=$outDirByCWD
fi
mkdir -p $outDir



echo "Running HOMER Peak Finder...: $chip"

# STEP 3: Run HOMER Peak Finder
if [ "$input" = "none" ]
then
	echo -e "findPeaks $chip -o $outDir/${name}.homer.txt $optStr &> $outDir/$name.homer.log" >&2
	findPeaks $chip -o $outDir/${name}.homer.txt $optStr 2>&1 | tee $outDir/$name.homer.log
else
	echo -e "findPeaks $chip -o $outDir/${name}.homer.txt -i $input $optStr &> $outDir/$name.homer.log" >&2
	findPeaks $chip -o $outDir/${name}.homer.txt -i $input $optStr 2>&1 | tee $outDir/$name.homer.log
fi

if [ $isFactor -eq 1 ];then
#	tbp=`grep "maximum tags considered per bp =" $peakFile | gawk '{printf "%d", $(NF)}'`
#	if [ "$tbp" == "" ];then
#	       echo "Error: Can't find tbp information" >&2
#	       exit 1
#	fi
#	if [ $tbp -eq 1 ];then
#		tDir=$chip/../TSVu
#	elif [ $tbp -gt 1 ];then
#		tDir=$chip/../TSV_tbp${tbp}
#		if [ ! -f $tDir/tagInfo.txt ];then
#			echo "TSV_tbp${tbp} does not exist. Generating.." >&2
#			if [ ! -d $chip/../TSVrd ];then
#				echo "Error: $chip/../TSVrd dos not exist." >&2
#				exit 1
#			fi
#			mkdir -p $tDir
#			makeTagDirectory $tDir -tbp $tbp -d $chip/../TSVrd
#		fi
#	else
#		echo "Error: Invalid tbp information." >&2
#		exit 1
#	fi

	postHomerPeak.sh $outDir/${name}.homer.txt $chip
	nPeak=`wc -l $outDir/${name}.homer.200bp.1rpm.bed | cut -d " " -f 1`
	nPeak0=`grep '^chr' $outDir/${name}.homer.txt | wc -l`
	echo -e "Final number of peak calls(>=1rpm, 200bp) $name: $nPeak (total $nPeak0)" >&2
else
	grep -v '^#' $outDir/${name}.homer.txt \
		| gawk '{w=$4-$3;rpkm=$6*1000/w; printf "%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%s\n", $2,$3,$4,$1,$6,$5,rpkm,w}'\
		> $outDir/${name}.homer.bed
fi

echo -e "" >&2
