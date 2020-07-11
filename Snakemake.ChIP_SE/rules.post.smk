
def get_downsample_depth(wildcards):
	if "DownSample" in samples:
		return samples.DownSample[samples.Name == wildcards.sampleName].tolist()[0]
	else:
		return None

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
		depth = get_downsample_depth
	message:
		"Downsampling... [{wildcards.sampleName}]"
	shell:
		"""
		module load ChIPseq/1.0
		bamDownsample.sh -n {params.depth} {input} > {output}
		"""

def get_align_bam(wildcards):
	# return ordered [ctrl , target] list.
	if "DownSample" in samples and samples.DownSample[samples.Name == wildcards.sampleName].tolist()[0] > 0:
		return alignDir + "/{sampleName}/align.ds.bam"
	else:
		return alignDir + "/{sampleName}/align.bam"

rule make_tagdir:
	input:
		get_align_bam
	output:
		directory(sampleDir + "/{sampleName}/TSV1")
	params:
		desDir = sampleDir + "/{sampleName}",
		name = "{sampleName}"
	message:
		"Making Homer tag directory... [{wildcards.sampleName}]"
	shell:
		"""
		module load ChIPseq/1.0
		mypipe.makeTagDir.sh -o {params.desDir} -n {params.name} -t 1 -c {chrRegexTarget} {input}
		"""


'''
## Not being used; replaced with get_input_tagdir for handling NULL ctrl
def get_input_name(sampleName):
	ctrlName = samples.Ctrl[samples.Name == sampleName]
	ctrlName = ctrlName.tolist()[0]
	return ctrlName
'''

'''
def get_input_tagdir(sampleName):
	ctrlName = samples.Ctrl[samples.Name == sampleName]
	ctrlName = ctrlName.tolist()[0]
	if ctrlName == "NULL":
		return "NULL"
	else:
		return sampleDir + "/" + ctrlName + "/TSV1"
'''

def get_peakcall_input(sampleName):
	ctrlName = samples.Ctrl[samples.Name == sampleName]
	ctrlName = ctrlName.tolist()[0]
	assert( len(crlName) == 1 )
	if ctrlName.upper() == "NULL":
		return [ sampleDir + "/" + sampleName + "/TSV1" ]
	else:
		return [ sampleDir + "/" + ctrlName + "/TSV1", sampleDir + "/" + samplename + "/TSV1" ]


rule call_peak_factor:
	input:
		lambda wildcards: get_peakcall_input(wildcards.sampleName)
		#chip = sampleDir + "/{sampleName}/TSV1",
		#ctrl = lambda wildcards: get_input_tagdir(wildcards.sampleName)
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed"
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	params:
		desDir = sampleDir + "/{sampleName}",
		mask = peak_mask,
		optStr = lambda wildcards, input: "-i" if len(input)>1 else ""
	shell:
		"""
		module load ChIPseq/1.0
		chip.peakCallFactor.sh -o {params.desDir}/HomerPeak.factor -m {params.mask} -s "-size 200" {params.optStr} {input}
		"""
#		chip.peakCallFactor.sh -o {params.desDir}/HomerPeak.factor -i {input.ctrl} -m {params.mask} -s "-size 200" {input}

'''
rule call_peak_factor_no_ctrl:
	input:
		chip = sampleDir + "/{sampleName}/TSV1",
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor.noCtrl/peak.exBL.1rpm.bed"
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	params:
		desDir = sampleDir + "/{sampleName}",
		mask = peak_mask
	shell:
		"""
		module load ChIPseq/1.0
		chip.peakCallFactor.sh -o {params.desDir}/HomerPeak.factor.noCtrl -m {params.mask} -s "-size 200" {input}
		"""
'''

rule call_peak_histone:
	input:
		lambda wildcards: get_peakcall_input(wildcards.sampleName)
		#chip = sampleDir + "/{sampleName}/TSV1",
		#ctrl = lambda wildcards: get_input_tagdir(wildcards.sampleName)
	output:
		sampleDir + "/{sampleName}/HomerPeak.histone/peak.exBL.bed"
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	params:
		desDir = sampleDir + "/{sampleName}",
		mask = peak_mask,
		optStr = lambda wildcards, input: "-i" if len(input)>1 else ""
	shell:
		"""
		module load ChIPseq/1.0
		chip.peakCallHistone.sh -o {params.desDir}/HomerPeak.histone -m {params.mask} {params.optStr} {input}
		"""


rule run_homermotif:
	input:
		sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed"
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed.all.noBG/homerResults.html"
	message:
		"Running Homer motif search... [{wildcards.sampleName}]"
	shell:
		"""
		module load Motif/1.0
		runHomerMotif.sh -g {genome} -s 200 -p 4 -b /data/limlab/Resource/Homer.preparse -o {sampleDir}/{wildcards.sampleName}/HomerPeak.factor {input}
		"""

'''
rule run_homermotif_no_ctrl:
	input:
		sampleDir + "/{sampleName}/HomerPeak.factor.noCtrl/peak.exBL.1rpm.bed"
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor.noCtrl/peak.exBL.1rpm.bed.all.noBG/homerResults.html"
	message:
		"Running Homer motif search... [{wildcards.sampleName}]"
	shell:
		"""
		module load Motif/1.0
		runHomerMotif.sh -g {genome} -s 200 -p 4 -b /data/limlab/Resource/Homer.preparse -o {sampleDir}/{wildcards.sampleName}/HomerPeak.factor.noCtrl {input}
		"""
'''

## RPM-scaled bigWig
rule make_bigwig:
	input:
		sampleDir + "/{sampleName}/TSV1"
	output:
		sampleDir + "/{sampleName}/igv.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	shell:
		"""
		module load ChIPseq/1.0
		chip.tagDirToBigWig.sh -g {chrom_size} -o {output} {input}
		"""
#	params:
#		memory = "%dG" % (  cluster["make_bigwig"]["memory"]/1000 - 1 )

