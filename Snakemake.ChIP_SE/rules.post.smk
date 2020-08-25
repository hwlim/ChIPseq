## Default values to main compatibility with previously precessed data
if 'Homer_tbp' not in locals():
	Homer_tbp=1

if 'downsampleDir' not in locals():
	downsampleDir="1.2.Align.downSample"

if 'doDedup' not in locals():
	doDedup=False

if 'dedupDir' not in locals():
	dedupDir="1.3.Align.dedup"

if 'bigWigDir_avg' not in locals():
	bigWigDir_avg="3.1.bigWig_avg"



def get_downsample_depth(wildcards):
	if "DownSample" in samples:
		return samples.DownSample[samples.Name == wildcards.sampleName].tolist()[0]
	else:
		return None

rule downsample_bam:
	input:
		alignDir + "/{sampleName}/align.bam"
	output:
		downsampleDir + "/{sampleName}/align.bam"
	params:
		depth = get_downsample_depth
	message:
		"Downsampling... [{wildcards.sampleName}]"
	shell:
		"""
		module load ChIPseq/1.0
		bamDownsample.sh -n {params.depth} {input} > {output}
		"""


def get_align_bam_for_dedup(wildcards):
	# return ordered [ctrl , target] list.
	if "DownSample" in samples and samples.DownSample[samples.Name == wildcards.sampleName].tolist()[0] > 0:
		return downsampleDir + "/{sampleName}/align.bam"
	else:
		return alignDir + "/{sampleName}/align.bam"

rule dedup_align:
	input:
		get_align_bam_for_dedup
	output:
		dedupDir + "/{sampleName}/align.bam"
	message:
		"Deduplicating... [{wildcards.sampleName}]"
	params:
		memory = "%dG" % ( cluster["dedup_align"]["memory"]/1000 - 2 )
	shell:
		"""
		module load CnR/1.0
		cnr.dedupBam.sh -m {params.memory} -o {output} -r {input}
		"""



def get_align_bam_for_tagdir(wildcards):
	# return ordered [ctrl , target] list.
	if doDedup:
		srcDir = dedupDir
	else:
		if "DownSample" in samples and samples.DownSample[samples.Name == wildcards.sampleName].tolist()[0] > 0:
			srcDir = downsampleDir
		else:
			srcDir = alignDir	
	return srcDir + "/{sampleName}/align.bam"

rule make_tagdir:
	input:
		get_align_bam_for_tagdir
	output:
		directory(sampleDir + "/{sampleName}/TSV")
	params:
		name = "{sampleName}"
	message:
		"Making Homer tag directory... [{wildcards.sampleName}]"
	shell:
		"""
		module load ChIPseq/1.0
		ngs.alignToTagDir.sh -o {output} -t {Homer_tbp} -c {chrRegexTarget} {input}
		drawAutoCorrplot.r -t {params.name} -o {output}/Autocorrelation.png {output}
		echo "{params.name}" > {sampleDir}/{wildcards.sampleName}/info.txt
		"""
#		mypipe.makeTagDir.sh -o {params.desDir} -n {params.name} -t {Homer_tbp} -c {chrRegexTarget} {input}


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
		return sampleDir + "/" + ctrlName + "/TSV"
'''

def get_input_name(sampleName):
	ctrlName = samples.Ctrl[samples.Name == sampleName]
	assert( len(ctrlName) == 1 )
	ctrlName = ctrlName.tolist()[0]
	return ctrlName

def get_peakcall_input(sampleName):
	ctrlName = get_input_name(sampleName)
	if ctrlName.upper() == "NULL":
		return [ sampleDir + "/" + sampleName + "/TSV" ]
	else:
		return [ sampleDir + "/" + ctrlName + "/TSV", sampleDir + "/" + sampleName + "/TSV" ]

def get_peakcall_opt(sampleName):
	if "PeakOpt" not in samples:
		return ""
	else:
		optStr = samples.PeakOpt[samples.Name == sampleName]
		assert( len(optStr) == 1 )
		optStr = optStr.tolist()[0]
		if optStr == "NULL":
			return ""
		else:
			return optStr

## NOTE: "-tbp 0" is implicitly set within chip.peakCallHistone.sh 
rule call_peak_factor:
	input:
		lambda wildcards: get_peakcall_input(wildcards.sampleName)
		#chip = sampleDir + "/{sampleName}/TSV",
		#ctrl = lambda wildcards: get_input_tagdir(wildcards.sampleName)
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed"
	message:
		"Calling TF peaks... [{wildcards.sampleName}]"
	params:
		desDir = sampleDir + "/{sampleName}",
		mask = peak_mask,
		optStr = lambda wildcards, input:( "-s \"" + get_peakcall_opt(wildcards.sampleName) + "\"" + " -i" ) if len(input)>1 else "-s \"" + get_peakcall_opt(wildcards.sampleName) + "\""
	shell:
		"""
		module load ChIPseq/1.0
		chip.peakCallFactor.sh -o {params.desDir}/HomerPeak.factor/peak -m {params.mask} -s "-size 200" {params.optStr} {input}
		"""
#		chip.peakCallFactor.sh -o {params.desDir}/HomerPeak.factor -i {input.ctrl} -m {params.mask} -s "-size 200" {input}


## NOTE: "-tbp 0" is implicitly set within chip.peakCallHistone.sh 
rule call_peak_histone:
	input:
		lambda wildcards: get_peakcall_input(wildcards.sampleName)
		#chip = sampleDir + "/{sampleName}/TSV",
		#ctrl = lambda wildcards: get_input_tagdir(wildcards.sampleName)
	output:
		sampleDir + "/{sampleName}/HomerPeak.histone/peak.exBL.bed"
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	params:
		desDir = sampleDir + "/{sampleName}",
		mask = peak_mask,
		optStr = lambda wildcards, input:( "-s \"" + get_peakcall_opt(wildcards.sampleName) + "\"" + " -i" ) if len(input)>1 else "-s \"" + get_peakcall_opt(wildcards.sampleName) + "\""
	shell:
		"""
		module load ChIPseq/1.0
		chip.peakCallHistone.sh -o {params.desDir}/HomerPeak.histone/peak -m {params.mask} {params.optStr} {input}
		"""


rule run_homer_motif:
	input:
		sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed"
	output:
		sampleDir + "/{sampleName}/Motif/Homer.all/homerResults.html"
	message:
		"Running Homer motif search... [{wildcards.sampleName}]"
	shell:
		"""
		module load Motif/1.0
		runHomerMotifSingle.sh -g {genome} -s 200 -p 4 -b /data/limlab/Resource/Homer.preparse \
			-o {sampleDir}/{wildcards.sampleName}/Motif/Homer.all {input}
		"""


rule run_meme_motif_rand5k:
	input:
		sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed"
	output:
		sampleDir + "/{sampleName}/Motif/MEME.random5k/meme-chip.html"
	message:
		"Running MEME-ChIP motif search for random 5k TSS peaks [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load MotifMEME/1.0
		runMemeChipSingle.sh -g {genomeFa} -s 200 -p 4 -r 5000 -d ~/bin/Motif/MEME_DB/Merged_By_Lim.meme \
			-o {sampleDir}/{wildcards.sampleName}/Motif/MEME.random5k {input}
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
		sampleDir + "/{sampleName}/TSV"
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

##### Average bigWig files
def get_bigwig_rep(groupName, srcDir):
	repL = samples.Name[samples.Group == groupName].tolist()
	return map(lambda x: srcDir + "/" + x + "/igv.bw", repL)

rule make_bigwig_avg:
	input:
		lambda wildcards: get_bigwig_rep(wildcards.groupName, sampleDir)
	output:
		bigWigDir_avg + "/{groupName}.avg.bw"
	message:
		"Making average bigWig files... [{wildcards.groupName}]"
	params:
		memory = "9G"
	shell:
		"""
		module load ChIPseq/1.0
		N=`ls {input} | wc -l`
		if [ $N -gt 1 ];then
			makeBigWigAverage.sh -g {chrom_size} -m {params.memory} -o {output} {input}
		else
			cp -v {input} {output}
		fi
		"""

def get_input_group(groupName):
	sampleName = samples.Name[samples.Group == groupName].tolist()
	assert( len(sampleName) >= 1 )
	sampleName = sampleName[0]

	ctrlName = samples.Ctrl[samples.Name == sampleName].tolist()
	assert( len(ctrlName) >= 1 )
	ctrlName = ctrlName[0]

	ctrlGroup = samples.Group[samples.Name == ctrlName].tolist()
	assert( len(ctrlGroup) == 1 )
	ctrlGroup = ctrlGroup[0]
	return ctrlGroup

rule make_bigwig_avg_subinput:
	input:
		chip = bigWigDir_avg + "/{groupName}.avg.bw",
		ctrl = lambda wildcards: bigWigDir_avg + "/" + get_input_group(wildcards.groupName) + ".avg.bw",
	output:
		bigWigDir_avg + "/{groupName}.avg.subInput.bw"
	message:
		"Making input-subtracted average bigWig files... [{wildcards.groupName}]"
	params:
		memory = "9G"
	#	memory = "%dG" % (  cluster["make_bigwig_subtract"]["memory"]/1000 - 1 )
	shell:
		"""
		module load CnR/1.0
		bigWigSubtract.sh -g {chrom_size} -m {params.memory} -t -1000 {output} {input.chip} {input.ctrl}
		"""
