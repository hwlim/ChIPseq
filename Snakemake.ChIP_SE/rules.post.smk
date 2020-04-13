
rule make_tagdir:
	input:
		alignDir+"/{sampleName}/align.bam"
	output:
		directory("{sampleName}/TSV1")
	params:
		name = "{sampleName}"
	message:
		"Making Homer tag directory... [{wildcards.sampleName}]"
	shell:
		"""
		module load ChIPseq/1.0
		mypipe.makeTagDir.sh -o {params.name} -n {params.name} -t 1 -c {chrRegexTarget} {input}
		"""

def get_input_name(sampleName):
	# return ordered [ctrl , target] list.
	ctrlName = samples.Ctrl[samples.Name == sampleName]
	ctrlName = ctrlName.tolist()[0]
	return ctrlName

rule call_peak_factor:
	input:
		chip="{sampleName}/TSV1",
		ctrl=lambda wildcards: get_input_name(wildcards.sampleName) + "/TSV1"
	output:
		"{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed"
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	params:
		name = "{sampleName}",
		mask = peak_mask
	shell:
		"""
		module load ChIPseq/1.0
		chip.peakCallFactor.dev.sh -o {params.name}/HomerPeak.factor -i {input.ctrl} -m {params.mask} -s "-size 200" {input}
		"""

rule call_peak_histone:
	input:
		chip="{sampleName}/TSV1",
		ctrl=lambda wildcards: get_input_name(wildcards.sampleName) + "/TSV1"
	output:
		"{sampleName}/HomerPeak.histone/peak.exBL.bed"
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	params:
		name = "{sampleName}",
		mask = peak_mask
	shell:
		"""
		module load ChIPseq/1.0
		chip.peakCallHistone.dev.sh -o {params.name}/HomerPeak.histone -i {input.ctrl} -m {params.mask} {input}
		"""


## RPM-scaled bigWig
rule make_bigwig:
	input:
		chip="{sampleName}/TSV1"
	output:
		"{sampleName}/igv.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	params:
		memory = "%dG" % (  cluster["make_bigwig"]["memory"]/1000 - 1 )
	shell:
		"""
		module load CnR/1.0
		chip.tagDirToBigWig.dev.sh -g {chrom_size} -m {params.memory} -o {output} {input}
		"""

