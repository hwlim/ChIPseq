#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh

if [ $# -lt 1 ];then
	srcConfig=./config.pipeline.PE.txt
else
	srcConfig=$1
fi
assertFileExist $srcConfig

echo -e "Loading data configuration file:\n\tSource = $srcConfig" >&2
configL=( ` cat $srcConfig | grep -v ^# | cut -f 4,5,6,7,8` )
# FastQ1/ FastQ2 / Genome / AlignOption / AlignmentPrefix

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



#optDefault="-q -S --best --strata --fullref -m 1 -p $thread --chunkmbs 512"
optDefault="-q --un /dev/null -p $thread"

for (( i=0;i<${#configL[@]};i=$i+5 ))
do
	fqFile1=${configL[$i]}
	fqFile2=${configL[$i+1]}
	genome=${configL[$i+2]}
	optStr=`echo ${configL[$i+3]} | sed 's/NULL//g' | sed 's/:/ /g'`
	desFile=${configL[$i+4]}
	name=${desFile%.*}

	optStr="${optDefault} ${optStr}"

	fq1=${readDir}/$fqFile1
	fq2=${readDir}/$fqFile2
	des=${desDir}/$desFile
	echo -e "Running Bowtie Paired-End" >&2
	echo -e "Source Files" >&2
	echo -e "  1: $fq1" >&2
	echo -e "  2: $fq2" >&2
	assertFileExist $fq1 $fq2

	echo -e "Parameter Setting" >&2
	echo -e "\tGenome = $genome" >&2
	echo -e "\toptStr = \"$optStr\"" >&2
	echo -e "\tdes = ${des}" >&2


	log=${desDir}/${name}.log

	echo -e "1) Running Bowtie" >&2
#	echo -e "bowtie $optStr $genome -1 $fq1 -2 $fq2" >&2
#	( bowtie $optStr $genome -1 $fq1 -2 $fq2 \
	echo -e "bowtie2 $optStr -x ${BT2_HOME}/index/$genome -1 $fq1 -2 $fq2" >&2
	( bowtie2 $optStr -x ${BT2_HOME}/index/$genome -1 $fq1 -2 $fq2 \
		| samtools view -bS - > ${des} ) 2>&1 | tee $log

#	echo -e "2) Sorting ${desUnsorted}" >&2
#	samtools sort -m 2G -T ${desDir}/${name} -O bam ${desUnsorted} > ${desSorted}

#	echo -e "3) Indexing $desSorted" >&2
#	samtools index ${desSorted}

#	rm ${desUnsorted}
#	exit 1
done
