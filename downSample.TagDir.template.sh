#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh

configL=(
mlds_500_4_ACAGTG/TSVrd	15M	15000000
mlds_500_4_CAGATC/TSVrd	15M	15000000
mlds_500_4_CGATGT/TSVrd	15M	15000000
mlds_500_4_CTTGTA/TSVrd	15M	15000000
mlds_500_4_GCCAAT/TSVrd	15M	15000000
mlds_500_4_TGACCA/TSVrd	15M	15000000
mlds_500_5_ACAGTG/TSVrd	15M	15000000
mlds_500_5_CAGATC/TSVrd	15M	15000000
mlds_500_5_CGATGT/TSVrd	15M	15000000
mlds_500_5_CTTGTA/TSVrd	15M	15000000
mlds_500_5_GCCAAT/TSVrd	15M	15000000
mlds_500_5_TGACCA/TSVrd	15M	15000000
mlds_500_6_ACAGTG/TSVrd	15M	15000000
mlds_500_6_CAGATC/TSVrd	15M	15000000
mlds_500_6_CGATGT/TSVrd	15M	15000000
mlds_500_6_CTTGTA/TSVrd	15M	15000000
mlds_500_6_GCCAAT/TSVrd	15M	15000000
mlds_500_6_TGACCA/TSVrd	15M	15000000
mlds_500_7_ACAGTG/TSVrd	15M	15000000
mlds_500_7_CAGATC/TSVrd	15M	15000000
mlds_500_7_CGATGT/TSVrd	15M	15000000
mlds_500_7_CTTGTA/TSVrd	15M	15000000
mlds_500_7_GCCAAT/TSVrd	15M	15000000
mlds_500_7_TGACCA/TSVrd	15M	15000000
mlds_605_6_ACAGTG/TSVrd	15M	15000000
mlds_605_6_CAGATC/TSVrd	15M	15000000
)

#Pol2
configL=(
mlds_629_8_ACAGTG/TSVrd 20M     20000000
mlds_629_8_CAGATC/TSVrd 20M     20000000
mlds_629_8_CGATGT/TSVrd 20M     20000000
mlds_629_8_CTTGTA/TSVrd 20M     20000000
mlds_629_8_GCCAAT/TSVrd 20M     20000000
mlds_629_8_TGACCA/TSVrd 20M     20000000
)

#Pol2
configL=(
mlds_629_8_ACAGTG/TSVrd	15M	15000000
mlds_629_8_CAGATC/TSVrd	15M	15000000
mlds_629_8_CGATGT/TSVrd	15M	15000000
mlds_629_8_CTTGTA/TSVrd	15M	15000000
mlds_629_8_GCCAAT/TSVrd	15M	15000000
mlds_629_8_TGACCA/TSVrd	15M	15000000
mlds_646_3_ACAGTG/TSVrd	15M	15000000
mlds_646_3_CAGATC/TSVrd	15M	15000000
mlds_646_3_CGATGT/TSVrd	15M	15000000
mlds_646_3_CTTGTA/TSVrd	15M	15000000
mlds_646_3_GCCAAT/TSVrd	15M	15000000
mlds_646_3_TGACCA/TSVrd	15M	15000000
)
#mkdir -p "$desDir"

for (( i=0;i<${#configL[@]};i=$i+3 ))
do
	srcTagDir=${configL[$i]}
	desSuffix=${configL[$i+1]}
	sampleSize=${configL[$i+2]}

	dataDir=${srcTagDir%\/TSVrd}

	TSVrd=${dataDir}/TSVrd.${desSuffix}
	TSVu=${dataDir}/TSVu.${desSuffix}

	tempBed=${dataDir}/__temp__.bed

	isDirExist ${srcTagDir}
	echo -e "Sampling:\n\tSrc: $srcTagDir\n\tDes: $TSVrd / $TSVu\n\tSize: $sampleSize" >&2
#	continue
	tagDir2bed.pl $srcTagDir -separate | shuf -n $sampleSize > ${tempBed}

	mkdir -p $TSVrd
	makeTagDirectory ${TSVrd} ${tempBed} -format bed 2>&1 | tee ${TSVrd}/TSVrd.log

	mkdir -p $TSVu
	makeTagDirectory ${TSVu} -d ${TSVrd} -tbp 1 2>&1 | tee ${TSVu}/TSVu.log

	rm ${tempBed}
done

