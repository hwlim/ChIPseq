########$$$$$$$###################################
### Snakemake rules for CUT&RUN Data Processing
### for Processing pooled replicates
###
### - except for "rule all"
### 
###		Written by Hee Woong Lim
########$$$$$$$###################################


def get_frag_replicate_names(groupName):
	repL = sampleAll.Name[sampleAll.Group == groupName].tolist()
	assert len(repL) > 0, "No replicates found for %s" % groupName
	return repL

## fragment pooling: frag.ctr
rule pool_replicate_frag_ctr:
	input:
		lambda wildcards: map(lambda x: frag_ctr_rep + "/" + x + ".frag.bed.gz", get_frag_replicate_names(wildcards.groupName))
	output:
		frag_ctr_pool + "/{groupName}.frag.bed.gz"
	params:
		memory= "%dG" % ( cluster["pool_replicate_frag_ctr"]["memory"]/1000 - 2 )
	message:
		"Pooling replicates... [{wildcards.groupName}]"
	shell:
		"""
		module load CnR/1.0
		zcat {input} | sort -S {params.memory} -k1,1 -k2,2n -k3,3n | gzip > {output}
		"""


rule pool_replicate_frag:
	input:
		lambda wildcards: map(lambda x: frag_ctr_rep + "/" + x + ".frag.bed.gz", get_frag_replicate_names(wildcards.groupName))
	output:
		frag_pool + "/{groupName}.frag.bed.gz"
	params:
		memory= "%dG" % ( cluster["pool_replicate_frag"]["memory"]/1000 - 2 )
	message:
		"Pooling replicates... [{wildcards.groupName}]"
	shell:
		"""
		module load CnR/1.0
		zcat {input} | sort -S {params.memory} -k1,1 -k2,2n -k3,3n | gzip > {output}
		"""
