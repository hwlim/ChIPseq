#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh

configL=( ` cat config.pipeline.txt | grep -v ^# | cut -f 7,8,9` )
# Alignment.RAW / Alignment.Norm / Alignment.Format


srcDir=1.Align.RAW
desDir=2.Align.Norm
mkdir -p "$desDir"

for (( i=0;i<${#configL[@]};i=$i+3 ))
do
	src=${srcDir}/${configL[$i]}
	des=${desDir}/${configL[$i+1]}
	field=${configL[$i+2]}
	assertFileExist ${src}

	echo -e "Reformatting input alignment file" >&2
	echo -e "\tSrc = $src" >&2
	echo -e "\tDes = $des" >&2
	echo -e "\tFormat = $field" >&2
	normalizeAlign.sh -f ${field} -sc ${src} > ${des}
	echo -e "" >&2
done

