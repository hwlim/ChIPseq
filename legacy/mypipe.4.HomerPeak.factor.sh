#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT


srcConfig=config.HomerPeak.factor.txt
assertFileExist $srcConfig

while IFS='' read -r line || [[ -n "$line" ]];
do
	[[ "$line" =~ ^#.*$ ]] &&  continue
	[[ "$line" = "" ]] && continue
	N_field=`echo "$line" | gawk -F "\t" '{ print NF }'`
	if [ "$N_field" -ne 6 ];then
		echo -e "Error: Invalid number of columns, $N_field" >&2
		exit 1
	fi  

	#############################################
	# Parameters 
	chip=`echo "$line" | cut -d $'\t' -f 1`
	input=`echo "$line" | cut -d $'\t' -f 2`
	name=`echo "$line" | cut -d $'\t' -f 3`
	desDir=`echo "$line" | cut -d $'\t' -f 4`
	optStr=`echo "$line" | cut -d $'\t' -f 5 | sed 's/,/ /g'`
	genome=`echo "$line" | cut -d $'\t' -f 6`

	blacklist=~/Research/Common_Data/${genome}/ENCODE-blacklist.bed
	#assertFileExist $blacklist
	isDirExist $chip

	#############################################
	# Output files
	des=${desDir}/${name}.homer.txt
	desBed0=${desDir}/${name}.homer.bed
	desBed200=${desDir}/${name}.homer.200bp.bed
	desBed200_exBL=${desDir}/${name}.homer.200bp.exBL.bed
	#desBed200_1rpm=${desDir}/${name}.homer.200bp.1rpm.bed
	desBed200_exBL_1rpm=${desDir}/${name}.homer.200bp.exBL.1rpm.bed
	log=${desDir}/${name}.homer.log

	#############################################
	# Peak calling
	echo -e "##########################################" >&2
	echo -e "Running Homer findPeaks: Factor mode" >&2
	echo -e "  Name = $name" >&2
	echo -e "  ChIP = $chip" >&2
	echo -e "  Input = $input" >&2
	echo -e "  desDir = $desDir" >&2
	echo -e "  optStr = $optStr" >&2
	echo -e "" >&2

	mkdir -p $desDir
	if [ "$optStr" = NULL ];then
		optStr=""
	fi

	if [ "$input" = NULL ];then
		echo -e "  1) findPeaks $chip -o ${desDir}/${name}.homer.txt $optStr -style factor" >&2
		findPeaks $chip -o ${des} $optStr -style factor 2>&1 | tee ${log}
	else
		isDirExist $input
		echo -e "  1) findPeaks $chip -o ${des} -i $input $optStr -style factor" >&2
		findPeaks $chip -o ${des} -i $input $optStr -style factor 2>&1 | tee ${log}
	fi
	pos2bed.pl $des | grep -v ^# | grep -v ^$ > $desBed0

	#############################################
	# Peak post-processing
	extendBed.sh -sq -w 200 $desBed0 > $desBed200
	echo -e "  2) Normalized tag counting (RPM)" >&2
	countTagsBed2.sh -q -ra $desBed200 $chip \
		| gawk '{ printf "%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,$8,$6 }' \
		| sort -k5,5nr \
		> __temp__.$$.bed
	mv __temp__.$$.bed $desBed200

	echo -e "  3) Filtering ENCODE blacklist regions" >&2
	if [ -f $blacklist ];then
		intersectBed -a $desBed200 -b $blacklist -v > $desBed200_exBL
	else
		echo -e "No black list file exists, pass" >&2
		cp $desBed200 $desBed200_exBL
	fi

	cat $desBed200_exBL | gawk '$5 > 1' > $desBed200_exBL_1rpm
	N=`cat $desBed200_exBL_1rpm | wc -l`
	echo -e "  >> Final number of peaks (>1RPM/exBL): $N" >&2
	echo -e "" >&2

done < $srcConfig

