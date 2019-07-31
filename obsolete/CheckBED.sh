#!/usr/bin/env bash

export LC_ALL=C
dataPath='./BED'

echo -e "FileName\tTotalTagCnt\tForward\tReverse\tF/R" > checkBED.info
for data in $(ls $dataPath/*.bed); do
	filename=${data##*\/}
	echo "Checking $filename"
	filename=${filename%%.bed}
	result=`gawk 'BEGIN {FR=0; RV=0}{if($6=="+")FR++; else RV++;} END{OFS="\t"; print FR+RV,FR,RV,FR/RV}' $data`
	echo -e "$filename\t$result" >> checkBED.info
done
unset LC_ALL
