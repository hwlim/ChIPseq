#!/usr/bin/env bash

nJob=30
totalWaitTime="48:00"
timestamp=$(date +%Y%m%d_%H%M%S)

if [ ! -f diag.pdf ];then
	module load python3/3.6.3
	snakemake --dag | dot -Tpdf > diag.pdf
fi

#module load python3/3.6.3
#snakemake -np
#exit 0
mkdir -p logs
mkdir -p logs.bsub
bsub -W ${totalWaitTime} -eo logs.bsub/bsub.${timestamp}.err -oo logs.bsub/bsub.${timestamp}.out \
	"module load python3/3.6.3
	snakemake -j $nJob \
		--latency-wait 60 \
		--cluster-config cluster.yml \
		--cluster 'bsub -W {cluster.walltime} -n {cluster.cpu} -M {cluster.memory} -J $$.{cluster.name} -R {cluster.resource} -eo {cluster.error} -oo {cluster.output}'"

