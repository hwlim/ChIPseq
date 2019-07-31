#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh

configL=(
#LN229_H3K27ac_Glu10mM_pool/TSV.fl150	"kwjl_0969_5_GTCCGC/TSVu.fl150 kwjl_0969_5_GTGAAA/TSVu.fl150"
#LN229_H3K27ac_Glu1mM_Ac5mM_pool/TSV.fl150	"kwjl_0969_5_CCGTCC/TSVu.fl150 kwjl_0969_6_ATGTCA/TSVu.fl150"
#LN229_H3K27ac_Glu1mM_pool/TSV.fl150	"kwjl_0969_5_AGTCAA/TSVu.fl150 kwjl_0969_5_AGTTCC/TSVu.fl150"

LN229_H3K27ac_Glu10mM_pool/TSV.fl150.15M	"kwjl_0969_5_GTCCGC/TSVu.fl150.15M kwjl_0969_5_GTGAAA/TSVu.fl150.15M"
LN229_H3K27ac_Glu1mM_Ac5mM_pool/TSV.fl150.15M	"kwjl_0969_5_CCGTCC/TSVu.fl150.15M kwjl_0969_6_ATGTCA/TSVu.fl150.15M"
LN229_H3K27ac_Glu1mM_pool/TSV.fl150.15M	"kwjl_0969_5_AGTCAA/TSVu.fl150.15M kwjl_0969_5_AGTTCC/TSVu.fl150.15M"
)

configL=(
LN229_H3K4me1_Glu10mM_pool/TSV.fl150.15M	"kwjl_0969_6_ACTGAT/TSVu.fl150.15M kwjl_0969_6_ATTCCT/TSVu.fl150.15M"
LN229_H3K4me1_Glu1mM_Ac5mM_pool/TSV.fl150.15M	"kwjl_0969_5_CGTACG/TSVu.fl150.15M kwjl_0969_6_GAGTGG/TSVu.fl150.15M"
LN229_H3K4me1_Glu1mM_pool/TSV.fl150.15M	"kwjl_0969_6_GTGGCC/TSVu.fl150.15M kwjl_0969_6_TGACCA/TSVu.fl150.15M"

LN229_H3K4me1_Glu10mM_pool/TSV.fl150	"kwjl_0969_6_ACTGAT/TSVu.fl150 kwjl_0969_6_ATTCCT/TSVu.fl150"
LN229_H3K4me1_Glu1mM_Ac5mM_pool/TSV.fl150	"kwjl_0969_5_CGTACG/TSVu.fl150 kwjl_0969_6_GAGTGG/TSVu.fl150"
LN229_H3K4me1_Glu1mM_pool/TSV.fl150	"kwjl_0969_6_GTGGCC/TSVu.fl150 kwjl_0969_6_TGACCA/TSVu.fl150"
)



for (( i=0;i<${#configL[@]};i=$i+2 ))
do
	desDir=${configL[$i]}
	tDirL=( ${configL[$i+1]} )
	dataDir=`dirname $desDir`
	name=$dataDir

	isDirExist ${tDirL[@]}

	echo -e "Pooling: ${desDir}" >&2
	echo ${tDirL[@]} | gawk '{ for(i=1;i<=NF;i++) printf "\t%s\n", $i }' >&2

#	continue
	mkdir -p ${desDir}
	echo -e "$name" > ${dataDir}/info.txt
	makeTagDirectory ${desDir} -d ${tDirL[@]} -fragLength 150 2>&1 | tee ${desDir}/TSV.log
done

