#/usr/bin/env bash

# Convert bin or nrm directory into a wig file
# Usage:
#	BIN2WIG.sh [srcDir] [desFile] [binsize]

if [ $# -lt 3 ];then
	echo -e "Usage: BIN2WIG.sh [srcDir] [desFile] [binsize]"
	exit
fi

src=$1
des=$2
binsize=$3

if [ -f $des ]; then
	rm -i $des
fi

for file in `ls $src/chr*.*`
do
#	echo $file
	chr=${file%%.*}
	chr=${chr##*/}
	printf "fixedStep chrom=%s start=1 step=%d\n" $chr $binsize >> $des
	gz=${file##*.}
#	echo $chr $gz
	if [ $gz="gz" ]; then
		zcat $file >> $des
	else
		cat $file >> $des
	fi
done
