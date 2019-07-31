#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh

source ./config.ChIP.star.sh

## data list file
if [ $# -gt 0 ];then
	srcDataList=$1
	assertFileExist $1
else
	echo -e "Error: Data list file must be specified" >&2
	exit 1
fi

srcBase=${fqBase}
desBase=${alignDir}
isDirExist $srcBase
mkdir -p $desBase

echo -e "============================================"
echo -e "Running STAR" >&2
echo -e "SrcBase: $srcBase" >&2
echo -e "DesBase: $desBase" >&2
echo -e "Index: $starIndex" >&2
echo -e "============================================"
echo -e "" >&2
isDirExist $starIndex


if [ -z ${bamSort+x} ]; then
	bamSort=TRUE
fi
if [ "${bamSort}" = "TRUE" ];then
	sortOpt="--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000"
else
	sortOpt="--outSAMtype BAM Unsorted"
fi


while IFS='' read -r line || [[ -n "$line" ]];
do
	[[ "$line" =~ ^#.*$ ]] &&  continue
	[[ "$line" = "" ]] && continue
	N_field=`echo "$line" | gawk -F "\t" '{ print NF }'`
	if [ "$N_field" -lt 5 ];then
		echo -e "Error: Invalid number of columns: Must be >= 5, currently $N_field" >&2
		exit 1
	fi
	
	id=`echo "$line" | cut -d "	" -f 1`
	dataDir=`echo "$line" | cut -d "	" -f 2`
	name=`echo "$line" | cut -d "	" -f 3`
	fq1=`echo "$line" | cut -d "	" -f 4`
	fq2=`echo "$line" | cut -d "	" -f 5`

	desPrefix=${desBase}/${id}.
	desBAM=${desPrefix}bam

	
	fq1=${srcBase}/${fq1}
	assertFileExist $fq1
	if [ "${fq2}" != "NULL" ];then
		fq2=${srcBase}/${fq2}
		assertFileExist $fq2
	fi
	

	echo -e "############################################" >&2
	echo -e "Running STAR" >&2
	echo -e "  id = $id" >&2
	echo -e "  name = $name" >&2
	echo -e "  fq1 = ${fq1}" >&2
	echo -e "  fq2 = ${fq2}" >&2;
	echo -e "  desBAM = $desBAM" >&2


	if [ -f ${desBAM} ];then
		echo -e "Warning: ${desBAM} already exists, skip" >&2
		continue
	fi


	if [ "${fq2}" == "NULL" ];then
		echo -e "STAR --runMode alignReads
	$optStr_star
	--genomeDir ${starIndex}
	--readFilesIn <( zcat $fq1 )
	$sortOpt
	--outReadsUnmapped $outReadsUnmapped
	--genomeLoad LoadAndKeep
	--outFileNamePrefix ${desPrefix}" >&2
		STAR --runMode alignReads \
			$optStr_star \
			--genomeDir ${starIndex} \
			--readFilesIn <( zcat $fq1 ) \
			$sortOpt \
			--outReadsUnmapped $outReadsUnmapped \
			--genomeLoad LoadAndKeep \
			--outFileNamePrefix ${desPrefix}
	else
		echo -e "STAR --runMode alignReads
	$optStr_star
	--genomeDir ${starIndex}
	--readFilesIn <( zcat $fq1 ) <( zcat $fq2 )
	$sortOpt
	--outReadsUnmapped $outReadsUnmapped
	--genomeLoad LoadAndKeep
	--outFileNamePrefix ${desPrefix}" >&2
		STAR --runMode alignReads \
			$optStr_star \
			--genomeDir ${starIndex} \
			--readFilesIn <( zcat $fq1 ) <( zcat $fq2 ) \
			${sortOpt} \
			--outReadsUnmapped $outReadsUnmapped \
			--genomeLoad LoadAndKeep \
			--outFileNamePrefix ${desPrefix}
	fi

	if [ "$bamSort" = "TRUE" ];then
		mv ${desPrefix}Aligned.sortedByCoord.out.bam ${desBAM}
		echo -e "Indexing" >&2
		samtools index ${desBAM}
	else
		mv ${desPrefix}Aligned.out.bam ${desBAM}
	fi


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

	echo -e "" >&2
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

