#!/usr/bin/env bash

## Written by Hee Woong Lim
##
## Script to submit snakemake job to CCHMC lsf system

source $COMMON_LIB_BASE/commonBash.sh

nJob=50
totalWaitTime="48:00"
#timestamp=$(date +%Y%m%d_%H%M%S)
if [ -z ${CHIP_SEQ_PATH+x} ]; then
	echo -e "Error: Environment variable CHIP_SEQ_PATH is not defined" >&2
	exit 1
fi

config=${CHIP_SEQ_PATH}/Snakemake.SE/cluster.yml

assertFileExist $config
assertFileExist ./Snakefile



# if [ ! -f diag.pdf ];then
# 	module load python3/3.6.3
# 	module load graphviz/2.40.1
# 	snakemake --dag | dot -Tpdf > diag.pdf
# fi

#module load python3/3.6.3
#snakemake -np
#exit 0
mkdir -p logs
bsub -W ${totalWaitTime} -eo submit.err -oo submit.out -q rhel9 <<- EOF
	module purge
	module load anaconda3
	source activate snakemake-7.18.2
	export XDG_CACHE_HOME=/scratch/$USER/snakemake-cache
	snakemake -j $nJob \
		--latency-wait 60 \
		--cluster-config $config \
		--cluster 'bsub -q rhel9 -W {cluster.walltime} -n {cluster.cpu} -M {cluster.memory} -J $$.{cluster.name} -R {cluster.resource} -eo {cluster.error} -oo {cluster.output}'
	unset XDG_CACHE_HOME
EOF

