########$$$$$$$###################################
### Snakemake rules for CUT&RUN Data Processing
### for Processing pooled replicates
###
### - except for "rule all"
### 
###		Written by Hee Woong Lim
########$$$$$$$###################################

def get_frag_replicate(groupName, repDir):
	repL = sampleAll.Name[sampleAll.Group == groupName].tolist()
	return map(lambda x: repDir + "/" + x + ".frag.bed.gz", repL)

## fragment pooling: frag.ctr
rule pool_replicate_frag_ctr:
	input:
		get_frag_replicate("{groupName}", frag_ctr_rep)
	output:
		frag_ctr_pool + "/{groupName}.frag.bed.gz"
	params:
		memory= "%dG" % ( cluster["pool_replicate_frag"]["memory"]/1000 - 1 )
	message:
		"Pooling replicates... [{wildcards.groupName}]"
	shell:
		"""
		module load CnR/1.0
		zcat {input} | sort -S {params.memory} -k1,1 -k2,2n -k3,3n | gzip > {output}
		"""


rule pool_replicate_frag:
	input:
		get_frag_replicate("{groupName}", frag_rep)
	output:
		frag_pool + "/{groupName}.frag.bed.gz"
	params:
		memory= "%dG" % ( cluster["pool_replicate_frag"]["memory"]/1000 - 1 )
	message:
		"Pooling replicates... [{wildcards.groupName}]"
	shell:
		"""
		module load CnR/1.0
		zcat {input} | sort -S {params.memory} -k1,1 -k2,2n -k3,3n | gzip > {output}
		"""
