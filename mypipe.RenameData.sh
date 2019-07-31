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

for (( i=0;i<${#configL[@]};i=$i+5 ))
do
	dataDir=${configL[$i]}
	name=${configL[$i+1]}
	genome=${configL[$i+2]}
	srcFile=${configL[$i+3]}
	flag=${configL[$i+4]}

	genome=${genome%_c}

	isDirExist $dataDir
	infoFile=${dataDir}/info.txt
	assertFileExist $infoFile

	## info.txt

	oldName=`cat ${infoFile}`
	echo -e "Updating data annotation\n\tOld = ${oldName}\n\tNew = ${name}\n" >&2
	echo $name > ${dataDir}/info.txt

	## HomerPeak files
	homerExist=`ls -f ${dataDir}/HomerPeak*/* > /dev/null 2>&1; if [ "$?" == "0" ]; then echo "TRUE"; else echo "FALSE"; fi`
	if [ "${homerExist}" = "FALSE" ];then continue; fi

	echo -e "Renaming HomerPeak file names" >&2
#	srcL=( `ls ${dataDir}/HomerPeak*/*` )
	srcDirL=( `ls -d ${dataDir}/HomerPeak*/` )
	for srcDir in ${srcDirL[@]}
	do
		echo -e "\tProcessing ${srcDir}" >&2
		srcL=( `ls ${srcDir}/${oldName}*` )
		for src in ${srcL[@]}
		do
			srcFile=`basename $src`
			srcDir=`dirname $src`
			
			if [[ ${srcFile} == ${oldName}* ]];then
				des=${srcDir}/${srcFile/${oldName}/${name}}
#				echo -e "\t${src} => ${des}" >&2
				mv ${src} ${des}
			fi
		done
	done
done

