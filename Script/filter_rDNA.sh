#!/usr/bin/env bash

# command-line application template

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT

function printUsage {
	echo "Usage: `basename $0` (options) [fq1] ..." >&2
	echo "Removes rDNA reads and creates new fastq file with a suffix, ex_rDNA.fq.gz" >&2
	echo -e "Options:" >&2
        echo -e "\t-g <rDNA index>: Path to Bowtie1 index for rDNA" >&2
        echo -e "\t-o <outDir>: Output directory, default=<same with src file>" >&2
        echo -e "\t-p <thread>: Thread for processing, default=4" >&2
#        echo -e "\t-f : flag option" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
genome=""
desDir=""
thread=4
while getopts ":g:o:p:" opt; do
	case $opt in
		g)
			genome=$OPTARG
			;;
		o)
			desDir=$OPTARG
			;;
		p)
			thread=$OPTARG
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

if [ "$genome" = "" ];then
	echo "Error: genome must be specified (-g)" >&2
	printUsage
	exit 1
fi

shift $((OPTIND-1))
if [ $# -eq 0 ];then
	printUsage
	exit 1
fi

assertFileExist $@
#assertFileExist $genome
srcL=( $@ )


###################################
## main code
echo -e "Filtering rDNA reads" >&2
echo -e "  rDNA index = $genome" >&2
echo -e "" >&2


if [ "$desDir" != "" ];then
	mkdir -p $desDir
fi

for src in ${srcL[@]}
do
	srcFile=`basename $src`
	if [ "$desDir" = "" ];then
		desDir=`dirname $src`
	fi
	prefix=${srcFile%.gz}
	prefix=${srcFile%.fastq}
	prefix=${srcFile%.fq}
	des=${desDir}/${prefix}.ex_rDNA.fq.gz
	log=${des%.fq.gz}.log

	echo -e "Processing $src" >&2
	echo -e "  => $des" >&2
	if [ -f $des ];then
		echo -e "  $des already exits, pass" >&2
		continue
	fi

	optStr="-q --fullref -p $thread --chunkmbs 512"
	( zcat $src \
		| bowtie $optStr $genome "-" --un __temp__.$$.fq \
		> /dev/null ) 2>&1 | tee $log
	gzip -v __temp__.$$.fq
	mv __temp__.$$.fq.gz $des
	echo -e "" >&2
done

exit 0
