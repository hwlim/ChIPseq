#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT


srcConfig=config.HomerPeak.histone.txt
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
	assertFileExist $blacklist
	isDirExist $chip

	#############################################
	# Output files
	des=${desDir}/${name}.homer.txt
	desBed=${desDir}/${name}.homer.bed
	desBed_exBL=${desDir}/${name}.homer.exBL.bed
	log=${desDir}/${name}.homer.log

	#############################################
	# Peak calling
	echo -e "##########################################" >&2
	echo -e "Running Homer findPeaks: Histone mode" >&2
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
		echo -e "  1) findPeaks $chip -o ${desDir}/${name}.homer.txt $optStr -style histone" >&2
		findPeaks $chip -o ${des} $optStr -style histone 2>&1 | tee ${log}
	else
		isDirExist $input
		echo -e "  1) findPeaks $chip -o ${des} -i $input $optStr -style histone" >&2
		findPeaks $chip -o ${des} -i $input $optStr -style histone 2>&1 | tee ${log}
	fi

	echo -e "  2) Normalized tag counting" >&2
	grep -v ^# $des | gawk '{ printf "%s\t%s\t%s\t%s\t%s\t%s\n", $2,$3,$4,$1,$6,$5 }' > __temp__.$$.bed
	countTagsBed2.sh -q -ka __temp__.$$.bed $chip \
		| gawk '{ printf "%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,$8,$6 }' \
		| sort -k5,5nr \
		> $desBed

	echo -e "  3) Filtering ENCODE blacklist regions" >&2
	intersectBed -a $desBed -b $blacklist -v > $desBed_exBL

	N=`cat $desBed_exBL | wc -l`
	echo -e "  >> Final number of peaks: $N" >&2

done < $srcConfig
