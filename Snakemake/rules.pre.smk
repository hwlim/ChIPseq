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
		#fq1 = lambda wildcards: trimDir + "/" + samples.Id[samples.Name == wildcards.sampleName] + "_1.trim.fq.gz",
		#fq2 = lambda wildcards: trimDir + "/" + samples.Id[samples.Name == wildcards.sampleName] + "_2.trim.fq.gz"
	output:
		alignDir+"/{sampleName}/align.bam"
	message:
		"Aligning... [{wildcards.sampleName}]"
	params:
		index=star_index,
		option=star_option
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
			{input}
		"""
#		STAR --runMode alignReads --genomeDir {params.index} \
#			--genomeLoad NoSharedMemory \
#			--readFilesIn <( zcat {input.fq1} ) <( zcat {input.fq2} ) \
#			--runThreadN {threads} \
#			{params.option} \
#			--outFileNamePrefix __temp__.$$ 2> {log}
#		mv __temp__.$$Aligned.out.bam {output}"

rule filter_align:
	input:
		alignDir+"/{sampleName}/align.bam"
	output:
		filteredDir + "/{sampleName}.filtered.bam"
	message:
		"Filtering... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		cnr.filterBam.sh  -o {output} -c "{chrRegexAll}" {input}
		"""

rule dedup_align:
	input:
		filteredDir + "/{sampleName}.filtered.bam"
	output:
		dedupDir + "/{sampleName}.dedup.bam"
	message:
		"Deduplicating... [{wildcards.sampleName}]"
	params:
		memory = "%dG" % ( cluster["dedup_align"]["memory"]/1000 - 1 )
	shell:
		"""
		module load CnR/1.0
		cnr.dedupBam.sh -m {params.memory} -o {output} -r {input}
		"""

rule check_baseFreq:
	input:
		filteredDir + "/{sampleName}.filtered.bam"
	output:
		read1 = baseFreqDir + "/{sampleName}.filtered.R1.freq.line.png",
		read2 = baseFreqDir + "/{sampleName}.filtered.R2.freq.line.png"
	message:
		"Checking baseFrequency... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		bamToBed.separate.sh -o {baseFreqDir} {input}
		checkBaseFreq.plot.sh -g {genomeFa} -c {chrRegexTarget} -o {baseFreqDir} {baseFreqDir}/{wildcards.sampleName}.filtered.R1.bed.gz
		checkBaseFreq.plot.sh -g {genomeFa} -c {chrRegexTarget} -o {baseFreqDir} {baseFreqDir}/{wildcards.sampleName}.filtered.R2.bed.gz
		"""

rule make_fragment:
	input:
		dedupDir + "/{sampleName}.dedup.bam" if doDedup else filteredDir + "/{sampleName}.filtered.bam"
	output:
		fragDir + "/{sampleName}.frag.bed.gz"
	params:
		memory = "%dG" % ( cluster["make_fragment"]["memory"]/1000 - 1 )
	message:
		"Making fragment bed files... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		bamToFragment.sh -o {output} -l 150 -s -m {params.memory} {input}
		"""
