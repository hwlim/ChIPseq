#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh


trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo "Usage: `basename $0` (options) [bam1] ..." >&2
	echo "Description: Separate & print alignments from Read1 and Read2" >&2
	echo "Ouptut: Two bed (gzipped) files in the current directory
	- <bamFile>.R1.bed.gz
	- <bamFile>.R2.bed.gz" >&2
	echo -e "Options:" >&2
	echo -e "\t-o : destination directory, default=." >&2
        echo -e "\t-r : 1 (read 1 only) / 2 (read 2 only) / both. default: both" >&2
        echo -e "\t-v : verbose mode. default: off" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
readSelect=both
verbose=FALSE
desDir="."
while getopts ":o:r:v" opt; do
	case $opt in
		o)
			desDir=$OPTARG
			;;
		r)
			readSelect=$OPTARG
			;;
		v)
			verbose=TRUE
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
assertFileExist ${srcL[@]}

mkdir -p $desDir


if [ "$readSelect" != "1" ] && [ "$readSelect" != "2" ] && [ "$readSelect" != "both" ];then
	echo -e "Error: Invalid selection of reads ($readSelect). must be 1 / 2 / both" >&2
	exit 1
fi


if [ "$readSelect" = "1" ] || [ "$readSelect" = "both" ];then
	r1=1
else
	r1=0
fi
	
if [ "$readSelect" = "2" ] || [ "$readSelect" = "both" ];then
	r2=1
else
	r2=0
fi

for src in ${srcL[@]}
do

	srcFile=`basename $src`
	prefix=${srcFile%.bam}

	des1=${desDir}/${prefix}.R1.bed.gz
	tmp1=${TMPDIR}/__temp__.$$.1.bed
	des2=${desDir}/${prefix}.R2.bed.gz
	tmp2=${TMPDIR}/__temp__.$$.2.bed
	

	if [ "$verbose" = "TRUE" ];then
		echo -e "Processing $src" >&2
		if [ $r1 -eq 1 ];then 
			echo -e "  - $des1" >&2
		fi
		if [ $r2 -eq 1 ];then 
			echo -e "  - $des2" >&2
		fi
	fi
	
	bamToBed -i $src \
		| gawk 'BEGIN{
				r1='$r1'
				r2='$r2'
				des1="'${tmp1}'"
				des2="'${tmp2}'"
				if(r1==1) printf "" > des1
				if(r2==1) printf "" > des2
			}
			{
				if(r1==1 && $4 ~ /\/1$/) print $0 >> des1
				if(r2==1 && $4 ~ /\/2$/) print $0 >> des2
			}'
	
	if [ $r1 -eq 1 ];then
		gzip ${tmp1}
		mv ${tmp1}.gz $des1
	fi
	
	if [ $r2 -eq 1 ];then
		gzip ${tmp2}
		mv ${tmp2}.gz $des2
	fi
done
