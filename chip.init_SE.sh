#!/usr/bin/env bash

if [ -z ${COMMON_LIB_BASE+x} ]; then
	echo -e "Error: LimLabBase is not perperly set up" >&2
	exit 1
fi
source $COMMON_LIB_BASE/commonBash.sh

if [ -z ${CHIP_PATH+x} ]; then
	echo -e "Error: Environment variable CHIP_PATH is not defined" >&2
	exit 1
fi

echo -e "Initializing ChIP-seq:ChIP_SE analysis" >&2


cp -i -v ${CHIP_PATH}/Snakemake.SE/sample.tsv .
cp -i -v ${CHIP_PATH}/Snakemake.SE/Snakefile .
