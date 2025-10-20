#!/usr/bin/env bash

#######
## Template script to convert old structure of analysis folder to new sample-by-sample structure


source $COMMON_LIB_BASE/commonBash.sh

sampleDir=2.Samples
mkdir -p $sampleDir



nameL=( `tail -n +2 sample.tsv | grep -v ^# | cut -f 1` )

moveRename(){
	local src=$1
	local des=$2
	local desDir=`dirname $des`

	if [ -f $src ];then
		if [ -f $des ];then
			echo -e "Warning $des already exists; skip" >&2
		else
			mkdir -p $desDir
			echo -e "Moving $src -> $des" >&2
			mv -v $src $des
		fi
	else
		echo -e "Warning $src does not exists; skip" >&2
	fi
	
}

for name in ${nameL[@]}
do
	echo -e "Processing $name" >&2

	src=1.4.Frag/fragLenHist/${name}.png
	des=${sampleDir}/${name}/QC/fragLen.png
	moveRename $src $des

	src=1.4.Frag/fragLenHist/${name}.txt
	des=${sampleDir}/${name}/QC/fragLen.txt
	moveRename $src $des

	src=1.4.Frag/spikeinCnt/${name}.spikeCnt.txt
	des=${sampleDir}/${name}/QC/spikeCnt.txt
	moveRename $src $des

	if [ 1 -eq 0 ];then
		src=1.4.Frag/${name}.frag.bed.gz
		des=${sampleDir}/${name}/fragment.bed.gz
		moveRename $src $des

		src=2.1.1.BigWig.ctr.RPM/${name}.ctr.rpm.bw
		des=${sampleDir}/${name}/igv.ctr.rpm.bw
		moveRename $src $des

		src=2.1.2.BigWig.ctr.RPM.subInput/${name}.ctr.rpm.subInput.bw
		des=${sampleDir}/${name}/igv.ctr.rpm.subInput.bw
		moveRename $src $des

		src=2.2.1.BigWig.frag.RPM/${name}.frag.rpm.bw
		des=${sampleDir}/${name}/igv.frag.rpm.bw
		moveRename $src $des

		src=2.2.2.BigWig.frag.RPM.subInput/${name}.frag.rpm.subInput.bw
		des=${sampleDir}/${name}/igv.frag.rpm.subInput.bw
		moveRename $src $des

		src=2.3.1.BigWig.ctr.RPSM/${name}.ctr.rpsm.bw
		des=${sampleDir}/${name}/igv.ctr.rpsm.bw
		moveRename $src $des

		src=2.3.2.BigWig.ctr.RPSM.subInput/${name}.ctr.rpsm.subInput.bw
		des=${sampleDir}/${name}/igv.ctr.rpsm.subInput.bw
		moveRename $src $des

		src=2.4.1.BigWig.frag.RPSM/${name}.frag.rpsm.bw
		des=${sampleDir}/${name}/igv.frag.rpsm.bw
		moveRename $src $des

		src=2.4.2.BigWig.frag.RPSM.subInput/${name}.frag.rpsm.subInput.bw
		des=${sampleDir}/${name}/igv.frag.rpsm.subInput.bw
		moveRename $src $des
	fi
done
