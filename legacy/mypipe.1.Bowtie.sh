#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh

if [ $# -lt 1 ];then
	srcConfig=./config.pipeline.txt
else
	srcConfig=$1
fi
assertFileExist $srcConfig

echo -e "Loading data configuration file:\n\tSource = $srcConfig" >&2
set -f
configL=( ` cat $srcConfig | grep -v ^# | cut -f 4,5,6,7` )
set +f
# FastQ / Genome / AlignOption / Alignment.RAW

if [ -f ./config.pipeline.sh ];then
	source ./config.pipeline.sh
	if [ -z "$readDir" ] || [ -z "$thread" ] || [ -z "$alignDir" ] || [ -z "$procDir" ];then
		echo -e "Error: At least one of srcDir / thread / desDir is not specified" >&2
		exit 1
	fi
else
	readDir=0.Reads
	alignDir=1.Align
	procDir=.
	thread=4
fi

srcDir=$readDir
desDir=$alignDir
mkdir -p $desDir



optDefault="-q -S --best --strata --fullref -m 1 -p $thread --chunkmbs 512"

for (( i=0;i<${#configL[@]};i=$i+4 ))
do
	set -f
	fqFileL=( `echo ${configL[$i]} | sed 's/;/ /g'` )
	set +f
	genome=${configL[$i+1]}
	optStr=`echo ${configL[$i+2]} | sed 's/NULL//g' | sed 's/;/ /g'`
	desFile=${configL[$i+3]}
	name=${desFile%.bam}

	optStr="${optDefault} ${optStr}"

	srcL=""
	echo -e "=====================" >&2
	echo -e "Running Bowtie" >&2
	echo -e "=====================" >&2
	echo -e "Source Files" >&2
	for fqFile in ${fqFileL[@]}
	do
		src=${srcDir}/${fqFile}
		assertFileExist $src
		srcL="${srcL} $src"
	done
	for src in ${srcL[@]}
	do
		echo -e "\t$src" >&2
	done

	echo -e "Parameter Setting" >&2
	echo -e "\tGenome = $genome" >&2
	echo -e "\toptStr = \"$optStr\"" >&2
	echo -e "\tdes = ${desDir}/${desFile}" >&2


	log=${desDir}/${name}.log
	desUnsorted=${desDir}/${name}.unsorted.bam
	desSorted=${desDir}/${name}.bam

	assertFileExist $src
	echo -e "1) Running Bowtie" >&2
	( zcat $srcL \
		| bowtie $optStr $genome "-" \
		| samtools view -bS - > ${desUnsorted} ) 2>&1 | tee $log

	echo -e "2) Sorting ${desUnsorted}" >&2
	samtools sort -m 2G -T ${desDir}/${name} -O bam ${desUnsorted} > ${desSorted}

	echo -e "3) Indexing $desSorted" >&2
	samtools index ${desSorted}

	rm ${desUnsorted}
	echo -e "" >&2
#	exit 1
done
