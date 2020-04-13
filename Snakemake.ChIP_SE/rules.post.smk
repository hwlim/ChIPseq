
def get_downsample_depth(wildcards):
	return samples.DownSample[samples.Name == wildcards.sampleName].tolist()[0]

'''
def get_downsample_bam(wildcards):
	depth = samples.DownSample[samples.Name == wildcards.sampleName].tolist()[0]
	return alignDir + "/{sampleName}/align." + ( "ds%dM" % depth ) + ".bam"

def get_downsample_suffix(wildcards):
	depth = samples.DownSample[samples.Name == wildcards.sampleName].tolist()[0]
	return "ds%dM" % depth
'''

rule downsample_bam:
	input:
		alignDir + "/{sampleName}/align.bam"
	output:
		alignDir + "/{sampleName}/align.ds.bam"
	params:
		depth=get_downsample_depth
	message:
		"Downsampling... [{wildcards.sampleName}]"
	shell:
		"""
		module load ChIPseq/1.0
		bamDownsample.sh -n {params.depth} {input} > {output}
		"""

def get_align_bam(wildcards):
	# return ordered [ctrl , target] list.
	if "DownSample" in samples and samples.DownSample[samples.Name == wildcards.sampleName].tolist()[0] == 0:
		return alignDir+"/{sampleName}/align.bam"
	else:
		return alignDir+"/{sampleName}/align.ds.bam"

rule make_tagdir:
	input:
		get_align_bam
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
		"{sampleName}/TSV1"
	output:
		"{sampleName}/igv.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	shell:
		"""
		module load ChIPseq/1.0
		chip.tagDirToBigWig.dev.sh -g {chrom_size} -o {output} {input}
		"""
#	params:
#		memory = "%dG" % (  cluster["make_bigwig"]["memory"]/1000 - 1 )

