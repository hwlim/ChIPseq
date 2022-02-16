#!/usr/bin/env bash

if [ -z ${COMMON_LIB_BASE+x} ]; then
	echo -e "Error: LimLabBase is not perperly set up" >&2
	exit 1
fi
source $COMMON_LIB_BASE/commonBash.sh

CHIP_PATH=~/bin/Pipeline/Snakemake.PE
#if [ -z ${CUTLERY+x} ]; then
#	echo -e "Error: Environment variable CUTLERY is not defined" >&2
#	exit 1
#fi

echo -e "Initializing ChIP-seq:PE analysis" >&2


cp -i -v ${CHIP_PATH}/sample.tsv .
cp -i -v ${CHIP_PATH}/Snakefile .
