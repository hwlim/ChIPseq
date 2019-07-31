#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh

if [ $# -lt 1 ];then
	srcConfig=./config.pipeline.txt
else
	srcConfig=$1
fi
assertFileExist $srcConfig

echo -e "Loading data configuration file:\n\tSource = $srcConfig" >&2
configL=( ` cat $srcConfig | grep -v ^# | cut -f 2,3,5,7,8` )
# DataDir / Name / Genome / Alignment / PipelineFlag

if [ -f ./config.pipeline.sh ];then
	source ./config.pipeline.sh
	if [ -z "$readDir" ] || [ -z "$thread" ] || [ -z "$alignDir" ];then
		echo -e "Error: At least one of readDir / alignDir / procDir / thread  not specified" >&2
		exit 1
	fi
else
	readDir=0.Reads
	alignDir=1.Align
	procDir="."
	thread=4
fi

#srcDir=2.Align.Norm
srcDir=$alignDir
desBase=$procDir
for (( i=0;i<${#configL[@]};i=$i+5 ))
do
	dataDir=${procDir}/${configL[$i]}
	name=${configL[$i+1]}
	genome=${configL[$i+2]}
	genome=${genome%_c}
	srcFile=${configL[$i+3]}
	flag=${configL[$i+4]}

	src=${srcDir}/${srcFile}
	assertFileExist ${src}

	dataPipeline2.sh -n ${name} -g ${genome} -o ${dataDir} ${flag} ${src}
done

