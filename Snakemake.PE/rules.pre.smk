## default STAR module
if 'star_module' not in locals():
	star_module = "STAR/2.5"

#if 'splitDir' not in locals():
#	splitDir = "1.4.Align.split"



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
	if doTrim:
		return [trimDir + "/" + samples.Id[samples.Name == wildcards.sampleName].tolist()[0] + "_1.trim.fq.gz",
			trimDir + "/" + samples.Id[samples.Name == wildcards.sampleName].tolist()[0] + "_2.trim.fq.gz"]
	else:
		return [fastqDir + "/" + samples.Fq1[samples.Name == wildcards.sampleName].tolist()[0],
			fastqDir + "/" + samples.Fq2[samples.Name == wildcards.sampleName].tolist()[0]]
	
#	raise(ValueError("Unrecognized wildcard value for 'endedness': %s" % wildcards.endedness))

rule align_pe:
	input:
		get_fastq
	output:
		alignDir+"/{sampleName}/align.bam",
		alignDir+"/{sampleName}/align.bam.bai"
	message:
		"Aligning... [{wildcards.sampleName}]"
	params:
		index=star_index,
		option=star_option,
		star_module = star_module
	log:
		alignDir + "/{sampleName}/star.log"
	threads:
		cluster["align_pe"]["cpu"]
	shell:
		"""
		module load ChIPseq/1.0
		module load {params.star_module}
		star.align.sh -g {params.index} \
			-o {alignDir}/{wildcards.sampleName}/align \
			-t {threads} \
			-p '{params.option}' \
			-s \
			{input}
		"""
#		STAR --runMode alignReads --genomeDir {params.index} \
#			--genomeLoad NoSharedMemory \
#			--readFilesIn <( zcat {input.fq1} ) <( zcat {input.fq2} ) \
#			--runThreadN {threads} \
#			{params.option} \
#			--outFileNamePrefix __temp__.$$ 2> {log}
#		mv __temp__.$$Aligned.out.bam {output}"

#		alignDir=expand(alignDir+"/{sampleName}", sampleName=samples.Name.tolist()),

def get_align_dir(bamList):
	import os.path
	return list(map(lambda x: os.path.dirname(x), bamList ))

rule make_align_stat_table:
	input:
		expand(alignDir+"/{sampleName}/align.bam", sampleName=samples.Name.tolist())
	output:
		qcDir + "/alignStat.txt"
	params:
		inputDir = get_align_dir(expand(alignDir+"/{sampleName}/align.bam", sampleName=samples.Name.tolist()))
	message:
		"Creating alignment stat file"
	shell:
		"""
		module load ChIPseq/1.0
		star.getAlignStats.r {params.inputDir} > {output}
		"""

# ## Eliminate scaffold/random chromosomes
# ## Retaining regular chromosomes and spikein (optional) chromosomes
# rule filter_align:
# 	input:
# 		alignDir+"/{sampleName}/align.bam"
# 	output:
# 		bam = filteredDir + "/{sampleName}.filtered.bam",
# 		bai = filteredDir + "/{sampleName}.filtered.bam.bai"
# 	message:
# 		"Filtering... [{wildcards.sampleName}]"
# 	shell:
# 		"""
# 		module load Cutlery/1.0
# 		cnr.filterBam.sh  -o {output.bam} -c "{chrRegexAll}" {input}
# 		samtools index {output.bam}
# 		"""
