## default STAR module
if 'star_module' not in locals():
	star_module="STAR/2.5"


rule trim_pe:
	input:
		fq1 = lambda wildcards: fastqDir + "/" + samples.Fq1[samples.Id == wildcards.sampleId],
		fq2 = lambda wildcards: fastqDir + "/" + samples.Fq2[samples.Id == wildcards.sampleId]
	output:
		fq1 = trimDir + "/{sampleId}_1.trim.fq.gz",
		fq2 = trimDir + "/{sampleId}_2.trim.fq.gz"
	message:
		"Trimming... [{wildcards.sampleId}]"
	params:
		adapter = adapter,
		minLen = trim_minLen,
		minQual = trim_minQual
	log:
		trimDir + "/{sampleId}.trim.log"
	shell:
		"""
		# Note: Needs to be implemented as a quantum transaction
		cutadapt -a {params.adapter} -A {params.adapter} --minimum-length {params.minLen} -q {params.minQual} \
			-o __temp__.$$.1.fq.gz -p __temp__.$$.2.fq.gz {input.fq1} {input.fq2} 2>&1 | tee {log}
		mv __temp__.$$.1.fq.gz {output.fq1}
		mv __temp__.$$.2.fq.gz {output.fq2} 
		"""

def get_fastq(wildcards):
	#print(wildcards.sampleName)
	if samples.Fq2[samples.Name == wildcards.sampleName].tolist()[0].upper() == "NULL":
		if doTrim:
			return trimDir + "/" + samples.Id[samples.Name == wildcards.sampleName].tolist()[0] + "_1.trim.fq.gz"
		else:
			return fastqDir + "/" + samples.Fq1[samples.Name == wildcards.sampleName].tolist()[0]
	else:
		if doTrim:
			return [trimDir + "/" + samples.Id[samples.Name == wildcards.sampleName].tolist()[0] + "_1.trim.fq.gz",
				trimDir + "/" + samples.Id[samples.Name == wildcards.sampleName].tolist()[0] + "_2.trim.fq.gz"]
		else:
			return [fastqDir + "/" + samples.Fq1[samples.Name == wildcards.sampleName].tolist()[0],
				fastqDir + "/" + samples.Fq2[samples.Name == wildcards.sampleName].tolist()[0]]
	
#	raise(ValueError("Unrecognized wildcard value for 'endedness': %s" % wildcards.endedness))


rule align_star:
	input:
		get_fastq
		#fq1 = lambda wildcards: trimDir + "/" + samples.Id[samples.Name == wildcards.sampleName] + "_1.trim.fq.gz",
		#fq2 = lambda wildcards: trimDir + "/" + samples.Id[samples.Name == wildcards.sampleName] + "_2.trim.fq.gz"
	output:
		alignDir + "/{sampleName}/align.bam"
	message:
		"Aligning... [{wildcards.sampleName}]"
	params:
		index=star_index,
		option=star_option,
		star_module = star_module
	log:
		alignDir + "/{sampleName}/star.log"
	threads:
		cluster["align_star"]["cpu"]
	shell:
		"""
		module load ChIPseq/1.0
		module load {params.star_module}
		star.align.sh -g {params.index} \
			-o {alignDir}/{wildcards.sampleName}/align \
			-t {threads} \
			-p '{params.option}'
			-s \
			{input}
		"""

