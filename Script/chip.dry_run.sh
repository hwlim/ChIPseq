#!/usr/bin/env bash
## Script to perform snakemake dry-run

if [ -z ${CHIP_SEQ_PATH+x} ]; then
	echo -e "Error: Environment variable CHIP_SEQ_PATH is not defined. Refer to the initial setup section in the github page at https://github.com/hwlim/ChIPseq for instructions on setting up this enrironment variable." >&2
	exit 1
fi

module purge
module load anaconda3
source activate snakemake-7.18.2
export XDG_CACHE_HOME=/scratch/$USER/snakemake-cache

if [ -e "config.yml" ]; then
	echo -e "Performing dry-run in Cutlery default mode..." >&2
	snakemake -np -s ${CUTLERY}/Snakemake/Snakefile_Default
else
	echo -e "Performing dry-run in Cutlery advanced mode..." >&2
	snakemake -np 
fi

unset XDG_CACHE_HOME
conda deactivate
module purge
