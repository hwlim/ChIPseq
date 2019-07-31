#!/usr/bin/env bash


#suppressPackageStartupMessages(require('geneplotter', quiet=TRUE))
suppressPackageStartupMessages(require('RColorBrewer', quiet=TRUE))
#suppressPackageStartupMessages(require('MASS', quiet=TRUE))
#suppressPackageStartupMessages(require('robustbase', quiet=TRUE))
suppressPackageStartupMessages(require('optparse', quiet=TRUE))
#suppressPackageStartupMessages(require('KernSmooth', quiet=TRUE))


source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT

function printUsage {
	echo "****** Under development *******" >&2
	echo "Description: Perform phiX read filtering and make a new fastq file" >&2
	echo "Usage: `basename $0` (options) [fastq] ..." >&2
	echo -e "Options:" >&2
        echo -e "\t-o <desDir>: Destination directory, default=<same with the source file>" >&2
#        echo -e "\t-f : flag option" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi
exit 0


###################################
## option and input file handling
desDir=""
while getopts ":g:" opt; do
	case $opt in
		o)
			desDir=$OPTARG
			;;
		0)
			echo -e "Flat -0 set" >&2
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

#if [ "$genome" = "" ];then
#	echo "Error: genome must be specified (-g)" >&2
#	printUsage
#	exit 1
#fi

shift $((OPTIND-1))
if [ $# -eq 0 ];then
	printUsage
	exit 1
fi

assertFileExist $@



	echo -e "STAR --runMode alignReads
	$optStr_star
	--genomeDir ${starIndex}
	--readFilesIn ${fqL}
	--readFilesCommand zcat
	--outSAMstrandField intronMotif
	--outSAMtype BAM Unsorted
	--outFilterMultimapNmax ${MultimapNmax}
	--outReadsUnmapped ${outReadsUnmapped}
	--alignMatesGapMax 2000
	--genomeLoad LoadAndKeep
	--outFileNamePrefix ${desPrefix}" >&2
	STAR --runMode alignReads \
		$optStr_star \
		--genomeDir ${starIndex} \
		--readFilesIn ${fqL} \
		--readFilesCommand zcat \
		--outSAMstrandField intronMotif \
		--outSAMtype BAM Unsorted \
		--outFilterMultimapNmax ${MultimapNmax} \
		--outReadsUnmapped ${outReadsUnmapped} \
		--alignMatesGapMax 2000 \
		--genomeLoad LoadAndKeep \
		--outFileNamePrefix ${desPrefix}
	mv ${desPrefix}Aligned.out.bam ${desBAM}

#	echo -e "\nIndexing" >&2
#	n_mapped=`samtools view $desBAM | wc -l`
#	if [ $n_mapped -lt 1 ];then
#		echo -e "  Warning: No alignment reads, creating empty index"
#		touch ${desBAM}.bai
#	else
#		echo -e ">> samtools index ${desBAM}" >&2
#		samtools index ${desBAM}
#	fi
	echo -e "" >&2

	# unmapped reads
	if [ "${outReadsUnmapped}" = "Fastx" ];then
		echo -e "Organizing unmapped reads" >&2
		if [ -f ${desDir}/alignment.STARUnmapped.out.mate1 ];then
			gzip -v -c ${desDir}/alignment.STARUnmapped.out.mate1 \
				> ${desDir}/unmapped.mate1.fastq.gz
			rm ${desDir}/alignment.STARUnmapped.out.mate1
		fi
		if [ -f ${desDir}/alignment.STARUnmapped.out.mate2 ];then
			gzip -v -c ${desDir}/alignment.STARUnmapped.out.mate2 \
				> ${desDir}/unmapped.mate2.fastq.gz
			rm ${desDir}/alignment.STARUnmapped.out.mate2
		fi
		echo -e "" >&2
	fi
	
done < $srcDataList


echo -e "" >&2
set +o pipefail
indexBeingUsed=`ps aux | grep -w STAR | grep $starIndex | grep -v grep | wc -l`
if [ $indexBeingUsed -gt 0 ];then
	echo -e "STAR index ${starIndex} is being used by other job. Index will be kept in shared memory" >&2
else
	echo -e "STAR index ${starIndex} is not being used, will be removed from shared memory" >&2
	STAR --genomeDir ${starIndex} --genomeLoad Remove
fi

