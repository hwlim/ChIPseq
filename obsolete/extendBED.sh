#!/usr/bin/env bash
# extend bed regions to given length
for oldBED in $@ 
do
	newBED=${oldBED%%.bed.gz}ext.bed.gz
	echo "$oldBED --> $newBED"
	zcat $oldBED \
		| gawk '{OFS="\t";start=$2;end=$3; if($6=="+") end=start+33; else start=end-33; print $1,start,end,$4,$5,$6}' \
		| sort -k1,1 -k2,2n \
		| gzip -c > $newBED
done
