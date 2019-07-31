#!/usr/bin/env bash

for file in $(ls *.bed)
do
	echo "Processing $file"
	cat $file | grep -v track > temp.bed && mv -f temp.bed $file
done
