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
		module purge
		module load cutadapt/2.1.0
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
		alignDir+"/{sampleName}/align.sortByName.bam"
		#alignDir+"/{sampleName}/align.bam.bai"
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
		module purge
		module load ChIPseq/1.0
		module load {params.star_module}
		star.align.sh -g {params.index} \
			-o {alignDir}/{wildcards.sampleName}/align \
			-t {threads} \
			-p '{params.option}' \
			{input}
		mv {alignDir}/{wildcards.sampleName}/align.bam {alignDir}/{wildcards.sampleName}/align.sortByName.bam
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
		expand(alignDir+"/{sampleName}/align.sortByName.bam", sampleName=samples.Name.tolist())
	output:
		qcDir + "/alignStat.txt"
	params:
		inputDir = get_align_dir(expand(alignDir+"/{sampleName}/align.sortByName.bam", sampleName=samples.Name.tolist())),
		outPrefix = lambda wildcards, output: __import__("re").sub(".txt$","", output[0])
	message:
		"Creating alignment stat file"
	shell:
		"""
		module purge
		module load R/4.4.0
		star.getAlignStats.r -o {params.outPrefix} {params.inputDir}
		"""


rule csort_bam:
	input:
		bam=alignDir+"/{sampleName}/align.sortByName.bam"
	output:
		bam = alignDir+"/{sampleName}/align.bam",
		bai = alignDir+"/{sampleName}/align.bam.bai"
	message:
		"Sorting by coordinate... [{wildcards.sampleName}]"
	threads:
		cluster["csort_bam"]["cpu"]
	shell:
		"""
		module purge
		module load ChIPseq/1.0
		samtools sort -o {output.bam} -T ${{TMPDIR}}/csort_bam.{wildcards.sampleName}.${{RANDOM}} -@ {threads} -m 2G {input.bam}
		samtools index {output.bam}
		"""

##########################################
## Rule to handle multimapper using CSEM

## Run CSEM for postprocessing of multimapper
rule run_csem:
	input:
		bam = alignDir+"/{sampleName}/align.sortByName.bam"
	output:
		bam = alignDir+"/{sampleName}/CSEM/align.bam"
	message:
		"Running CSEM... [{wildcards.sampleName}]"
	params:
		#outDir = lambda wildcards, output: __import__("os").path.dirname(output[0]),
		prefix = lambda wildcards, output: __import__("re").sub(".bam$", "", output[0])
	threads:
		cluster["run_csem"]["cpu"]
	shell:
		"""
		module purge
		module load R/4.4.0
		module load bedtools/2.30.0
		module load csem_limlab/06272024
		ngs.run_csem.sh -o {params.prefix} -t {threads} -n {wildcards.sampleName} {input.bam}
		"""

## Unify CSEM bam file
rule unify_csem:
	input:
		alignDir+"/{sampleName}/CSEM/align.bam"
	output:
		alignDir+"/{sampleName}/CSEM/align.uniq.bam"
	message:
		"Unifying CSEM results... [{wildcards.sampleName}]"
	params:
		outDir = lambda wildcards, output: __import__("os").path.dirname(output[0]),
		prefix = lambda wildcards, output: __import__("re").sub(".bam$", "", output[0])
	shell:
		"""
		module purge
		module load ChIPseq/1.0
		ngs.unifyCSEM.py -o {output} {input}
		"""

# ## coordinate-sort unified CSEM bam file
# rule csort_csem:
# 	input:
# 		bam = alignDir+"/{sampleName}/CSEM/align.uniq.bam"
# 	output:
# 		bam = alignDir+"/{sampleName}/CSEM/align.uniq.sorted.bam",
# 		bai = alignDir+"/{sampleName}/CSEM/align.uniq.sorted.bam.bai"
# 	message:
# 		"Sorting by coordinate... [{wildcards.sampleName}]"
# 	threads:
# 		cluster["csort_csem"]["cpu"]
# 	shell:
# 		"""
# 		module load ChIPseq/1.0
# 		samtools sort -o {output.bam} -T ${{TMPDIR}}/csort_csem.{wildcards.sampleName}.${{RANDOM}} -@ {threads} -m 2G {input.bam}
# 		samtools index {output.bam}
# 		"""

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
# 		module load ChIPseq/1.0
# 		ngs.filterBam.sh  -o {output.bam} -c "{chrRegexAll}" {input}
# 		samtools index {output.bam}
# 		"""
