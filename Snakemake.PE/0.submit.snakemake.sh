#!/usr/bin/env bash

nJob=30
totalWaitTime="48:00"
timestamp=$(date +%Y%m%d_%H%M%S)
yml=~/bin/Pipeline/Snakemake/cluster.yml

if [ ! -f diag.pdf ];then
	module load python3/3.6.3
	module load graphviz/2.40.1
	snakemake --dag | dot -Tpdf > diag.pdf
fi

#module load python3/3.6.3
#snakemake -np
#exit 0
mkdir -p logs
bsub -W ${totalWaitTime} -eo bsub.${timestamp}.err -oo bsub.${timestamp}.out \
	"module load python3/3.6.3
	snakemake -j $nJob \
		--latency-wait 60 \
		--cluster-config $yml \
		--cluster 'bsub -W {cluster.walltime} -n {cluster.cpu} -M {cluster.memory} -J $$.{cluster.name} -R {cluster.resource} -eo {cluster.error} -oo {cluster.output}'"

