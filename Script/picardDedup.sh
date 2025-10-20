#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT


if [ $# -lt 1 ];then
	echo -e "Usage: picardMarkdup.sh <bamFile>" >&2
	echo -e "Description: Wrapper script for MarkDuplicates command in Picard tools" >&2
	echo -e "\tRemove duplicates and create new bam file with suffix .dedup.bam" >&2
	echo -e "\tAssumes input bam file is coordinate-sorted" >&2
	echo -e "Output: ${srcFile}.dedup.bam file is created in the current directory" >&2
	exit 1
fi

src=$1
assertFileExist $src
des=`basename $src`
des=${des%.bam}.dedup.bam
metric=${des%.bam}.metric

echo -e "Coordinate sorting by Picard" >&2
echo -e "  Src = $src" >&2
echo -e "  Des = $des" >&2
java -Xmx2000m -jar ${PICARD_HOME}/picard.jar MarkDuplicates \
	I=$src \
	O=__temp__.$$.bam \
	M=$metric \
	REMOVE_DUPLICATES=true 
mv __temp__.$$.bam $bam


