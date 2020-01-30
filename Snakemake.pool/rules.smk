########$$$$$$$###################################
### Snakemake rules for CUT&RUN Data Processing
### for Processing pooled replicates
###
### - except for "rule all"
### 
###		Written by Hee Woong Lim
########$$$$$$$###################################

## Replicate-pooling -> bam file
def get_frag_replicate(wildcards):
	repL = sampleAll.Name[sampleAll.Group == wildcards.groupName].tolist()
	return map(lambda x: fragDir_rep + "/" + x + ".frag.bed.gz", repL)

rule pool_replicate_frag:
	input:
		get_frag_replicate
	output:
		fragDir_pool + "/{groupName}.frag.bed.gz"
	params:
		memory= "%dG" % ( cluster["pool_replicate_frag"]["memory"]/1000 - 1 )
	message:
		"Pooling replicates... [{wildcards.groupName}]"
	shell:
		"""
		module load CnR/1.0
		zcat {input} | sort -S {params.memory} -k1,1 -k2,2n -k3,3n | gzip > {output}
		"""
