if 'Homer_tbp' not in locals():
	Homer_tbp=0

if 'fc_histone' not in locals():
	fc_histone=4

if 'spikein_chrom_size' not in locals():
	spikein_chrom_size="NULL"

########## Auxilary functions definition start #################

## Decide downsampling depth if any
def get_downsample_depth(wildcards):
	'''
	Get down sampling depth
	Default value is from 'downSampleN
	If 'DownSample' column exists in sample.tsv file, it has priority
	'''
	if "downSampleN" in locals():
		n = downSampleN
	else:
		n = 0

	if "DownSample" in samples:
		tmp = samples.DownSample[samples.Name == wildcards.sampleName].tolist()[0]
		if tmp > 0:
			n = tmp
	
	if isinstance(n, int) and n >= 0:
		return n
	else:
		raise RuntimeError('Invalid down sampling depth: %s' % n)

'''
## NOTE: These two functions cannot be used for Snakemake in normal situation because spikeCnt.txt file does not exist upfront in most of times
## Extract spikein % from spikein data file
def extract_spikeinfrac(src):
	# src="spikeCnt.txt"
	import pandas as pd
	import sys
	data = pd.read_csv(src, header=0, index_col="Name", sep="\t", comment="#", na_filter=False)
	assert("SpikeinFrac" in data.index)
	data = data.set_axis(["value"], axis=1)
	spikeinFrac = data["value"]["SpikeinFrac"]
	return spikeinFrac

## Calculate spikein ratio to multiply to fold-change criteria
def get_spikein_ratio(chip, ctrl):
	spikeChip = extract_spikeinfrac(sampleDir + "/" + chip + "/QC/spikeCnt.txt")
	spikeCtrl = extract_spikeinfrac(sampleDir + "/" + ctrl + "/QC/spikeCnt.txt")
	return spikeChip / spikeCtrl
'''

########## Auxilary functions definition end #################



## spikein information processing yielding
## - main read count from chromosomes starting with "chr"
## - spikein read count from chromosomes starting with spikein prefix
## - unknown read count from chomoromes starting with others
## - spikein % out of spike + main, not unknown
## - scaleFactor for RPSM normalization to multiply to raw read count
rule count_spikein:
	input:
		#fragDir + "/{sampleName}.frag.bed.gz"
		sampleDir + "/{sampleName}/fragment.bed.gz"
	output:
		sampleDir + "/{sampleName}/QC/spikeCnt.txt"
	message:
		"Counting spikein tags... [{wildcards.sampleName}]"
	shell:
		"""
		module load Cutlery/1.0
		ngs.countSpikein.sh -p {spikePrefix} -n {wildcards.sampleName} {input} > {output}
		"""

rule make_spikeintable:
	input:
		expand(sampleDir + "/{sampleName}/QC/spikeCnt.txt", sampleName=samples.Name.tolist())
	output:
		qcDir + "/spikein.txt"
	message:
		"Making spikein table..."
	params:
		outPrefix = lambda wildcards, output: __import__("re").sub(".txt$","", output[0])
	shell:
		"""
		module load Cutlery/1.0
		ngs.makeSpikeCntTable.r -o {params.outPrefix} {input}
		"""

## Draw a plot of fragment length distribution
## Automatically consider targer chromosomes only starting with "chr" excluding others such as "DM-chr"
rule get_fragLenHist:
	input:
		sampleDir + "/{sampleName}/fragment.bed.gz"
	output:
		txt = sampleDir + "/{sampleName}/QC/fragLen.txt",
		png = sampleDir + "/{sampleName}/QC/fragLen.png"
	message:
		"Checking fragment length... [{wildcards.sampleName}]"
	params:
		outPrefix = lambda wildcards, output: __import__("re").sub(".txt$","", output[0])
	shell:
		"""
		module load Cutlery/1.0
		ngs.fragLenHist.r -o {params.outPrefix} -n {wildcards.sampleName} {input}
		"""

## bigwig file: resized fragment in RPM scale
## Draw a plot of fragment length distribution
## Automatically consider targer chromosomes only starting with "chr" excluding others such as "DM-chr"
## and the same in other bigwig file rules
rule make_bigwig_ctr_rpm:
	input:
		sampleDir + "/{sampleName}/fragment.bed.gz"
	output:
		sampleDir + "/{sampleName}/igv.ctr.rpm.bw"
	message:
		"Making bigWig files, ctr.rpm ... [{wildcards.sampleName}]"
	params:
		memory = "%dG" % (  cluster["make_bigwig"]["memory"]/1000 - 2 )
	shell:
		"""
		module load Cutlery/1.0
		cnr.fragToBigWig.sh -g {chrom_size} -c "{chrRegexTarget}" -r 150 -m {params.memory} -o {output} {input}
		"""

## bigwig file: original fragment in RPM scale
rule make_bigwig_frag_rpm:
	input:
		sampleDir + "/{sampleName}/fragment.bed.gz"
	output:
		sampleDir + "/{sampleName}/igv.frag.rpm.bw",
	message:
		"Making bigWig files, frag.rpm ... [{wildcards.sampleName}]"
	params:
		memory = "%dG" % (  cluster["make_bigwig"]["memory"]/1000 - 2 )
	shell:
		"""
		module load Cutlery/1.0
		cnr.fragToBigWig.sh -g {chrom_size} -c "{chrRegexTarget}" -m {params.memory} -o {output} {input}
		"""


## bigwig file: resized fragment in RPSM scale
## scaling factor already calculated in spikeCnt.txt 
rule make_bigwig_ctr_rpsm:
	input:
		frag = sampleDir + "/{sampleName}/fragment.bed.gz",
		spikeinCnt = sampleDir + "/{sampleName}/QC/spikeCnt.txt"
	output:
		sampleDir + "/{sampleName}/igv.ctr.rpsm.bw",
	message:
		"Making bigWig files, ctr.rpsm ... [{wildcards.sampleName}]"
	params:
		memory = "%dG" % ( cluster["make_bigwig"]["memory"]/1000 - 2 )
	shell:
		"""
		module load Cutlery/1.0
		scaleFactor=`cat {input.spikeinCnt} | gawk '{{ if($1=="ScaleFactor") print $2 }}'`
		if [ "$scaleFactor" == "" ];then
			echo -e "Error: empty scale factor" >&2
			exit 1
		fi
		cnr.fragToBigWig.sh -g {chrom_size} -c "{chrRegexTarget}" -r 150 -m 5G -s $scaleFactor -o {output} {input.frag}
		"""

## bigwig file: original fragment in RPSM scale
rule make_bigwig_frag_rpsm:
	input:
		frag = sampleDir + "/{sampleName}/fragment.bed.gz",
		spikeinCnt = sampleDir + "/{sampleName}/QC/spikeCnt.txt"
	output:
		sampleDir + "/{sampleName}/igv.frag.rpsm.bw",
	message:
		"Making bigWig files, frag.rpsm ... [{wildcards.sampleName}]"
	params:
		memory = "%dG" % (  cluster["make_bigwig"]["memory"]/1000 - 2 )
	shell:
		"""
		module load Cutlery/1.0
		scaleFactor=`cat {input.spikeinCnt} | gawk '{{ if($1=="ScaleFactor") print $2 }}'`
		if [ "$scaleFactor" == "" ];then
			echo -e "Error: empty scale factor" >&2
			exit 1
		fi
		cnr.fragToBigWig.sh -g {chrom_size} -c "{chrRegexTarget}" -m 5G -s $scaleFactor -o {output} {input.frag}
		"""
#		scaleFactor=`cat {input.spikeinCnt} | gawk '$1=="'{wildcards.sampleName}'"' | gawk '{{ printf "%f", 100000/$3 }}'`


def get_ctrl_name(sampleName):
	# return ordered [ctrl , target] list.
	ctrlName = samples.Ctrl[samples.Name == sampleName]
	ctrlName = ctrlName.tolist()[0]
	""" Temporary debug code
	print("%s : %s" % (sampleName, ctrlName), file=sys.stderr)
	src=bigWigDir_ctr_RPM + "/" + ctrlName + ".ctr.rpm.bw"
	import os.path
	if os.path.isfile(src):
		print("%s exists" % ctrlName, file=sys.stderr)
	else:
		print("%s not exist" % ctrlName, file=sys.stderr)
	"""
	return ctrlName


rule make_bigwig_ctr_rpm_subinput:
	input:
		chip = sampleDir + "/{sampleName}/igv.ctr.rpm.bw",
		ctrl = lambda wildcards: sampleDir + "/" + get_ctrl_name(wildcards.sampleName) + "/igv.ctr.rpm.bw"
	output:
		sampleDir + "/{sampleName}/igv.ctr.rpm.subInput.bw",
	message:
		"Making bigWig files, ctr.rpm.subInput ... [{wildcards.sampleName}]"
	#params:
	#	memory = "%dG" % (  cluster["make_bigwig_subtract"]["memory"]/1000 - 1 )
	shell:
		"""
		module load Cutlery/1.0
		bigWigSubtract.sh -g {chrom_size} -m 5G -t -1000 {output} {input.chip} {input.ctrl}
		"""

rule make_bigwig_frag_rpm_subinput:
	input:
		chip = sampleDir + "/{sampleName}/igv.frag.rpm.bw",
		ctrl = lambda wildcards: sampleDir + "/" + get_ctrl_name(wildcards.sampleName) + "/igv.frag.rpm.bw"
	output:
		sampleDir + "/{sampleName}/igv.frag.rpm.subInput.bw",
	message:
		"Making bigWig files, frag.rpm.subInput ... [{wildcards.sampleName}]"
	#params:
	#	memory = "%dG" % (  cluster["make_bigwig_subtract"]["memory"]/1000 - 1 )
	shell:
		"""
		module load Cutlery/1.0
		bigWigSubtract.sh -g {chrom_size} -m 5G -t -1000 {output} {input.chip} {input.ctrl}
		"""


rule make_bigwig_ctr_rpsm_subinput:
	input:
		chip = sampleDir + "/{sampleName}/igv.ctr.rpsm.bw",
		ctrl = lambda wildcards: sampleDir + "/" + get_ctrl_name(wildcards.sampleName) + "/igv.ctr.rpsm.bw"
	output:
		sampleDir + "/{sampleName}/igv.ctr.rpsm.subInput.bw",
	message:
		"Making bigWig files, ctr.rpsm.subInput ... [{wildcards.sampleName}]"
	#params:
	#	memory = "%dG" % (  cluster["make_bigwig_subtract"]["memory"]/1000 - 1 )
	shell:
		"""
		module load Cutlery/1.0
		bigWigSubtract.sh -g {chrom_size} -m 5G -t -1000 {output} {input.chip} {input.ctrl}
		"""

rule make_bigwig_frag_rpsm_subinput:
	input:
		chip = sampleDir + "/{sampleName}/igv.frag.rpsm.bw",
		ctrl = lambda wildcards: sampleDir + "/" + get_ctrl_name(wildcards.sampleName) + "/igv.frag.rpsm.bw"
	output:
		sampleDir + "/{sampleName}/igv.frag.rpsm.subInput.bw",
	message:
		"Making bigWig files, frag.rpsm.subInput ... [{wildcards.sampleName}]"
	#params:
	#	memory = "%dG" % (  cluster["make_bigwig_subtract"]["memory"]/1000 - 1 )
	shell:
		"""
		module load Cutlery/1.0
		bigWigSubtract.sh -g {chrom_size} -m 5G -t -1000 {output} {input.chip} {input.ctrl}
		"""

'''
## RPM-scaled bigWig input-divided log2FC
rule make_bigwig_divide:
	input:
		get_bigwig_input
	output:
		bigWigDir_div + "/{sampleName}.divInput.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	#params:
	#	memory = "%dG" % (  cluster["make_bigwig_subtract"]["memory"]/1000 - 1 )
	shell:
		"""
		module load Cutlery/1.0
		bigWigDivide.sh -g {chrom_size} -m 5G -s log -a 1 -o {output} {input}
		"""

rule make_tagdir:
	input:
		fragDir + "/{sampleName}.frag.bed.gz"
	output:
		directory(homerDir + "/{sampleName}/TSV")
	params:
		name = "{sampleName}"
	message:
		"Making Homer tag directory... [{wildcards.sampleName}]"
	shell:
		"""
		module load Cutlery/1.0
		cnr.makeHomerDir.sh -c {chrRegexTarget} -o {output} -n {params.name} {input}
		"""


def get_hetchr_ctrl(sampleName):
	#print(sampleName)
	ctrlName = samples.Ctrl[samples.Name == sampleName]
	ctrlName = ctrlName.tolist()[0]
	return ctrlName
#	return [ fragDir + "/" + sampleName + ".frag.bed.gz", fragDir + "/" + ctrlName + ".frag.bed.gz" ]
'''

## TODO: The currrent version of Spikein does not precisely reflect spikein scaling factor
##      "Spikein" must be repladced with "SpikeinFrac"
## Need to double check if the "spikein" factor is correct (note that it is inverse of call_peak_hetchr_spikein_homer"
rule call_peak_hetchr:
	input:
		chip = sampleDir + "/{sampleName}/fragment.bed.gz",
		ctrl = lambda wildcards: sampleDir + "/" + get_ctrl_name(wildcards.sampleName) + "/fragment.bed.gz",
		spikeChip = sampleDir + "/{sampleName}/QC/spikeCnt.txt",
		spikeCtrl = lambda wildcards: sampleDir + "/" + get_ctrl_name(wildcards.sampleName) + "/QC/spikeCnt.txt"
	output:
		bed = peakDir + "/{sampleName}." + peakSuffix + ".bed",
		txt = peakDir + "/{sampleName}." + peakSuffix + ".txt.gz"
	params:
		outPrefix = lambda wildcards, output: __import__("re").sub(".bed$","", output[0])
	message:
		"Peak calling for heterochromatin by fold-change ... [{wildcards.sampleName}]"
	shell:
		"""
		module load ChIPseq
		spikeChip=`cat {input.spikeChip} | gawk '$1 == "Spikein"' | cut -f 2`
		spikeCtrl=`cat {input.spikeCtrl} | gawk '$1 == "Spikein"' | cut -f 2`
		if [ "$spikeChip" == "" ] || [ "$spikeCtrl" == "" ];then
			echo -e "Error: empty spikein factor" >&2
			exit 1
		fi
		spikein=`echo -e "${{spikeChip}}\t${{spikeCtrl}}" | gawk '{{ printf "%f", $2 / $1 }}'`
		hetChr.call.fc.sh -o {params.outPrefix} -g {chrom_size} -w {peakWindow} -s {peakStep} -a {peakAlpha} -k $spikein -t {peakFC} {input.chip} {input.ctrl}
		"""


rule call_peak_hetchr_spikein_homer:
	input:
		chip = sampleDir + "/{sampleName}/fragment.bed.gz",
		ctrl = lambda wildcards: sampleDir + "/" + get_ctrl_name(wildcards.sampleName) + "/fragment.bed.gz",
		spikeChip = sampleDir + "/{sampleName}/QC/spikeCnt.txt",
		spikeCtrl = lambda wildcards: sampleDir + "/" + get_ctrl_name(wildcards.sampleName) + "/QC/spikeCnt.txt"
	output:
		hetChr_homer_spike + "/{sampleName}.exBL.bed"
	params:
		outPrefix = lambda wildcards, output: __import__("re").sub(".exBL.bed$","", output[0])
	message:
		"Peak calling for heterochromatin by Homer/frag ... [{wildcards.sampleName}]"
	shell:
		"""
		module load ChIPseq
		spikeChip=`cat {input.spikeChip} | gawk '$1 == "Spikein"' | cut -f 2`
		spikeCtrl=`cat {input.spikeCtrl} | gawk '$1 == "Spikein"' | cut -f 2`
		if [ "$spikeChip" == "" ] || [ "$spikeCtrl" == "" ];then
			echo -e "Error: empty spikein factor" >&2
			exit 1
		fi
		spikein=`echo -e "${{spikeChip}}\t${{spikeCtrl}}" | gawk '{{ printf "%f", $1 / $2 }}'`
		hetChr.call.homer.sh -o {params.outPrefix} -c "{chrRegexTarget}" -s {peakWindow_homer} -d {minDist_homer} -k $spikein -f {peakFC_homer} -m {peak_mask} {input.chip} {input.ctrl}
		"""

## hetero chromatin peak calling using Homer
rule call_peak_hetchr_spikein_homer_ctr:
	input:
		chip = sampleDir + "/{sampleName}/fragment.bed.gz",
		ctrl = lambda wildcards: sampleDir + "/" + get_ctrl_name(wildcards.sampleName) + "/fragment.bed.gz",
		spikeChip = sampleDir + "/{sampleName}/QC/spikeCnt.txt",
		spikeCtrl = lambda wildcards: sampleDir + "/" + get_ctrl_name(wildcards.sampleName) + "/QC/spikeCnt.txt"
	output:
		hetChr_homer_spike_ctr + "/{sampleName}.exBL.bed",
	params:
		outPrefix = lambda wildcards, output: __import__("re").sub(".exBL.bed$","", output[0])
	message:
		"Peak calling for heterochromatin by Homer/ctr ... [{wildcards.sampleName}]"
	shell:
		"""
		module load ChIPseq
		spikeChip=`cat {input.spikeChip} | gawk '$1 == "Spikein"' | cut -f 2`
		spikeCtrl=`cat {input.spikeCtrl} | gawk '$1 == "Spikein"' | cut -f 2`
		if [ "$spikeChip" == "" ] || [ "$spikeCtrl" == "" ];then
			echo -e "Error: empty spikein factor" >&2
			exit 1
		fi
		spikein=`echo -e "${{spikeChip}}\t${{spikeCtrl}}" | gawk '{{ printf "%f", $1 / $2 }}'`
		hetChr.call.homer.sh -o {params.outPrefix} -c "{chrRegexTarget}" -r 150 -s {peakWindow_homer} -d {minDist_homer} -k $spikein -f {peakFC_homer} -m {peak_mask} -l 150 {input.chip} {input.ctrl}
		"""





################ Single-END style rules START #######################
def get_align_bam_for_tagdir(wildcards):
	# return ordered [ctrl , target] list.
	downDepth=get_downsample_depth(wildcards)
	if doDedup:
		srcDir = dedupDir
	else:
		#if "DownSample" in samples and samples.DownSample[samples.Name == wildcards.sampleName].tolist()[0] > 0:
		if downDepth > 0:
			srcDir = downsampleDir
		else:
			srcDir = alignDir	
	return srcDir + "/{sampleName}/align.bam"

rule make_tagdir_se:
	input:
		get_align_bam_for_tagdir
	output:
		directory(sampleDir + "/{sampleName}/TSV.SE")
	params:
		name = "{sampleName}",
		optStr = "" if "robustFragLen" not in locals() or robustFragLen == False else "-r"
	message:
		"Making Homer tag directory/SE... [{wildcards.sampleName}]"
	shell:
		"""
		module load ChIPseq/1.0
		ngs.alignToTagDir.sh -o {output} -t {Homer_tbp} -c {chrRegexTarget} {params.optStr} {input}
		drawHomerAutoCorr.r -t {params.name} -o {output}/Autocorrelation.png {output}
		echo "{params.name}" > {sampleDir}/{wildcards.sampleName}/info.txt
		"""


def get_peakcall_input_tagdir(sampleName):
	ctrlName = get_ctrl_name(sampleName)
	if ctrlName.upper() == "NULL":
		return [ sampleDir + "/" + sampleName + "/TSV.SE" ]
	else:
		return [ sampleDir + "/" + ctrlName + "/TSV.SE", sampleDir + "/" + sampleName + "/TSV.SE" ]

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

## NOTE: "-tbp 0" is implicitly set within chip.peakCallFactor.sh
## Peak size is fixed as 200bp
rule call_peak_factor:
	input:
		lambda wildcards: get_peakcall_input_tagdir(wildcards.sampleName)
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed"
	message:
		"Calling TF peaks/SE ... [{wildcards.sampleName}]"
	params:
		mask = peak_mask,
		outPrefix = lambda wildcards, output: __import__("re").sub(".exBL.1rpm.bed$","", output[0]),
		optStr = lambda wildcards, input:( "\"-size 200 " + get_peakcall_opt(wildcards.sampleName) + "\"" + " -i" ) if len(input)>1 else "\"-size 200 " + get_peakcall_opt(wildcards.sampleName) + "\""
	shell:
		"""
		module load ChIPseq/1.0
		chip.peakCallFactor.sh -o {params.outPrefix} -m {params.mask} -s {params.optStr} {input}
		"""
#		chip.peakCallFactor.sh -o {params.desDir}/HomerPeak.factor -i {input.ctrl} -m {params.mask} -s "-size 200" {input}


## NOTE: "-tbp 0" is implicitly set within chip.peakCallHistone.sh 
rule call_peak_histone:
	input:
		lambda wildcards: get_peakcall_input_tagdir(wildcards.sampleName)
	output:
		sampleDir + "/{sampleName}/HomerPeak.histone/peak.exBL.bed"
	message:
		"Calling histone peaks/SE ... [{wildcards.sampleName}]"
	params:
		mask = peak_mask,
		outPrefix = lambda wildcards, output: __import__("re").sub(".exBL.bed$","", output[0]),
		optStr = lambda wildcards, input:( "\"" + get_peakcall_opt(wildcards.sampleName) + "\"" + " -i" ) if len(input)>1 else "\"" + get_peakcall_opt(wildcards.sampleName) + "\""
	shell:
		"""
		module load ChIPseq/1.0
		chip.peakCallHistone.sh -o {params.outPrefix} -m {params.mask} -f {fc_histone} -s {params.optStr}  {input}
		"""


## homer histone peak calling considering spikein
rule call_peak_histone_spikein:
	input:
		tagDir = lambda wildcards: get_peakcall_input_tagdir(wildcards.sampleName),
		spikeChip = sampleDir + "/{sampleName}/QC/spikeCnt.txt",
		spikeCtrl = lambda wildcards: sampleDir + "/" + get_ctrl_name(wildcards.sampleName) + "/QC/spikeCnt.txt"
	output:
		sampleDir + "/{sampleName}/HomerPeak.histone.spikein/peak.exBL.bed"
	message:
		"Calling histone peaks/SE with spikein ... [{wildcards.sampleName}]"
	params:
		#spikeFactor = lambda wildcards: get_spikein_ratio(wildcards.sampleName, get_ctrl_name(wildcards.sampleName)),
		outPrefix = lambda wildcards, output: __import__("re").sub(".exBL.bed$","", output[0]),
		optStr = lambda wildcards, input:( "\"" + get_peakcall_opt(wildcards.sampleName) + "\"" + " -i" ) if len(input)>1 else "\"" + get_peakcall_opt(wildcards.sampleName) + "\""
	shell:
		"""
		module load ChIPseq/1.0
		spikeChip=`cat {input.spikeChip} | gawk '$1 == "SpikeinFrac"' | cut -f 2`
		spikeCtrl=`cat {input.spikeCtrl} | gawk '$1 == "SpikeinFrac"' | cut -f 2`
		if [ "$spikeChip" == "" ] || [ "$spikeCtrl" == "" ];then
			echo -e "Error: empty spikein factor" >&2
			exit 1
		fi
		spikeFactor=`echo -e "${{spikeChip}}\t${{spikeCtrl}}" | gawk '{{ printf "%f", $1 / $2 }}'`
		chip.peakCallHistone.sh -o {params.outPrefix} -m {peak_mask} -f {fc_histone}  -k $spikeFactor -s {params.optStr} {input.tagDir}
		"""

## homer TF peak calling considering spikein
rule call_peak_factor_spikein:
	input:
		tagDir = lambda wildcards: get_peakcall_input_tagdir(wildcards.sampleName),
		spikeChip = sampleDir + "/{sampleName}/QC/spikeCnt.txt",
		spikeCtrl = lambda wildcards: sampleDir + "/" + get_ctrl_name(wildcards.sampleName) + "/QC/spikeCnt.txt"
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor.spikein/peak.exBL.1rpm.bed"
	message:
		"Calling histone peaks/SE with spikein ... [{wildcards.sampleName}]"
	params:
		#spikeFactor = lambda wildcards: get_spikein_ratio(wildcards.sampleName, get_ctrl_name(wildcards.sampleName)),
		outPrefix = lambda wildcards, output: __import__("re").sub(".exBL.1rpm.bed$","", output[0]),
		optStr = lambda wildcards, input:( "\"" + get_peakcall_opt(wildcards.sampleName) + "\"" + " -i" ) if len(input)>1 else "\"" + get_peakcall_opt(wildcards.sampleName) + "\""
	shell:
		"""
		module load ChIPseq/1.0
		spikeChip=`cat {input.spikeChip} | gawk '$1 == "SpikeinFrac"' | cut -f 2`
		spikeCtrl=`cat {input.spikeCtrl} | gawk '$1 == "SpikeinFrac"' | cut -f 2`
		if [ "$spikeChip" == "" ] || [ "$spikeCtrl" == "" ];then
			echo -e "Error: empty spikein factor" >&2
			exit 1
		fi
		spikeFactor=`echo -e "${{spikeChip}}\t${{spikeCtrl}}" | gawk '{{ printf "%f", $1 / $2 }}'`
		chip.peakCallFactor.sh -o {params.outPrefix} -m {peak_mask} -f 4 -k $spikeFactor -s {params.optStr} {input}
		"""


rule run_homer_motif:
	input:
		sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed"
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor/Motif/Homer.all/homerResults.html"
	message:
		"Running Homer motif search... [{wildcards.sampleName}]"
	params:
		outPrefix =  lambda wildcards, output: __import__("os").path.dirname(output[0])
	shell:
		"""
		module load Motif/1.0
		runHomerMotifSingle.sh -g {genome} -s 200 -p 4 -b /data/limlab/Resource/Homer.preparse \
			-o {params.outPrefix} {input}
		"""
#		outPrefix=`dirname {output}`


rule run_meme_motif_rand5k:
	input:
		sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed"
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor/Motif/MEME.random5k/meme-chip.html"
	message:
		"Running MEME-ChIP motif search for random 5k TSS peaks [{wildcards.sampleName}]"
	params:
		outPrefix =  lambda wildcards, output: __import__("os").path.dirname(output[0])
	shell:
		"""
		module purge
		module load MotifMEME/1.0
		outPrefix=`dirname {output}`
		runMemeChipSingle.sh -g {genomeFa} -s 200 -p 4 -r 5000 -d ~/bin/Motif/MEME_DB/Merged_By_Lim.meme \
			-o {params.outPrefix} {input}
		"""


rule draw_peak_heatmap_factor:
	input:
		bed = sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed",
		bw = sampleDir + "/{sampleName}/igv.ctr.rpm.bw"
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor/heatmap.exBL.1rpm.png"
	message:
		"Drawing peak profile heatmap... [{wildcards.sampleName}]"
	params:
		outPrefix = lambda wildcards, output: __import__("re").sub(".png$","", output[0])
	shell:
		"""
		module load Cutlery/1.0
		drawBigWigHeatmap.r -t {wildcards.sampleName} -m 0,0.5,2,0.5 -w 2000 -b 20 -s 3,6 \
			-o {params.outPrefix} \
			{input.bed} {input.bw}
		"""


rule draw_peak_heatmap_histone:
	input:
		bed = sampleDir + "/{sampleName}/HomerPeak.histone/peak.exBL.bed",
		bw = sampleDir + "/{sampleName}/igv.ctr.rpm.bw"
	output:
		sampleDir + "/{sampleName}/HomerPeak.histone/heatmap.exBL.png"
	message:
		"Drawing peak profile heatmap... [{wildcards.sampleName}]"
	params:
		outPrefix = lambda wildcards, output: __import__("re").sub(".png$","", output[0])
	shell:
		"""
		module load Cutlery/1.0
		drawBigWigHeatmap.r -t {wildcards.sampleName} -m 0,0.5,2,0.5 -w 10000 -b 20 -s 3,6 \
			-o {params.outPrefix}\
			{input.bed} {input.bw}
		"""



##### Rules for spike-in bigwig files ########################

## bigwig file: resized fragment in RPM scale
## Draw a plot of fragment length distribution
## Automatically consider targer chromosomes only starting with "chr" excluding others such as "DM-chr"
## and the same in other bigwig file rules
rule make_bigwig_ctr_rpm_spike:
	input:
		sampleDir + "/{sampleName}/fragment.bed.gz"
	output:
		sampleDir + "/{sampleName}/spikein.igv.ctr.rpm.bw"
	message:
		"Making bigWig files, ctr.rpm ... [{wildcards.sampleName}]"
	params:
		memory = "%dG" % (  cluster["make_bigwig"]["memory"]/1000 - 2 )
	shell:
		"""
		module load Cutlery/1.0
		ngs.fragToSpikeBigWig.sh -g {spikein_chrom_size} -p "{spikePrefix}" -r 150 -m {params.memory} -o {output} {input}
		"""







################ Single-END style END #######################





'''
## Peak calling using Homer
## Rule from CUT&RUN pipeline

def get_peakcall_input(wildcards):
	# return ordered [ctrl , target] list. if no ctrl, simply [target].
	ctrlName = samples.Ctrl[samples.Name == wildcards.sampleName]
	ctrlName = ctrlName.tolist()[0]
	if ctrlName.upper() == "NULL":
		return [ homerDir + "/" + wildcards.sampleName + "/TSV" ]
	else:
		return [ homerDir + "/" + ctrlName + "/TSV", homerDir + "/" + wildcards.sampleName + "/TSV" ]

rule call_peaks:
	input:
		get_peakcall_input
		#target 	= homerDir + "/{sampleName}/TSV",
		#ctrl 	= "test"
		#ctrl 	= "get_ctrl"
	output:
		homerDir + "/{sampleName}/HomerPeak/peak.homer.exBL.1rpm.bed"
	params:
		mask = peak_mask,
		peakDir = homerDir + "/{sampleName}/HomerPeak",
		optStr = lambda wildcards, input: "-i" if len(input)>1 else ""
	message:
		"Peak calling using Homer... [{wildcards.sampleName}]"
	shell:
		"""
		module load Cutlery/1.0
		cnr.peakCallTF.sh -o {params.peakDir} -m {params.mask} -s \"-fragLength 100\" {params.optStr} {input}
		"""
'''





################################################################
## *** Note *****
## get_scalefactor function is deterministic, i.e Snakemake run this function to generate command lines to submit
## Therefore, when spikein.txt does not exist, this rule will invoke error.
## This rule must be revised the command take spikein.txt as an input directly not using the get_scalefactor function
#  
#####################################################
## Scaled BigWig by Spike-in using raw read counts (not RPM)
#def get_scalefactor(wildcards):
#	# return ordered [ctrl , target] list.
#	spikeinTable = pd.read_csv(spikeinCntDir + "/spikein.txt", sep="\t", comment="#", na_filter=False)
#	if not spikeinTable.Sample.is_unique:
#		print( "Error: Sample column spikein.txt is not unique")
#		sys.exit()
#	return [ spikeinTable.ScaleFactor[spikeinTable.Sample == wildcards.sampleName].tolist()[0] ]

## Raw read count scale + spike-in scaled
'''
rule make_bigwig_scaled:
	input:
		bed = fragDir + "/{sampleName}.frag.bed.gz",
		spikeinCnt = spikeinCntDir + "/spikein.txt"
	output:
		bigWigScaledDir + "/{sampleName}.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
#	params:
#		memory = "%dG" % ( cluster["make_bigwig"]["memory"]/1000 - 1 ),
#		scaleFactor = get_scalefactor
	shell:
		"""
		module load Cutlery/1.0
		scaleFactor=`cat {input.spikeinCnt} | gawk '$1=="'{wildcards.sampleName}'"' | cut -f 6`
		if [ "$scaleFactor" == "" ];then
			echo -e "Error: empty scale factor" >&2
			exit 1
		fi
		cnr.fragToBigWig.sh -g {chrom_size} -m 5G -s $scaleFactor -o {output} {input.bed}
		"""
#		cnr.fragToBigWig.sh -g {chrom_size} -m {params.memory} -s {params.scaleFactor} -o {output} {input.bed}

def get_bigwig_scaled_input(wildcards):
	# return ordered [ctrl , target] list.
	ctrlName = samples.Ctrl[samples.Name == wildcards.sampleName]
	ctrlName = ctrlName.tolist()[0]
	return [ bigWigScaledDir + "/" + wildcards.sampleName + ".bw", bigWigScaledDir + "/" + ctrlName + ".bw" ]

rule make_bigwig_scaled_subtract:
	input:
		get_bigwig_scaled_input
	output:
		bigWigScaledDir_sub + "/{sampleName}.scaled.subInput.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	#params:
	#	memory = "%dG" % (  cluster["make_bigwig_subtract"]["memory"]/1000 - 1 )
	shell:
		"""
		module load Cutlery/1.0
		bigWigSubtract.sh -g {chrom_size} -m 5G -t -1000 {output} {input}
		"""

rule make_bigwig_scaled_divide:
	input:
		get_bigwig_scaled_input
	output:
		bigWigScaledDir_div + "/{sampleName}.scaled.divInput.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	#params:
	#	memory = "%dG" % (  cluster["make_bigwig_subtract"]["memory"]/1000 - 1 )
	shell:
		"""
		module load Cutlery/1.0
		bigWigDivide.sh -g {chrom_size} -m 5G -s log -a 1 -o {output} {input}
		"""
'''


'''
##### Average bigWig files
def get_bigwig_rep_sub(wildcards):
	repL = samples.Name[samples.Group == wildcards.groupName].tolist()
	return map(lambda x: bigWigDir_sub + "/" + x + ".subInput.bw", repL)

rule make_bigwig_sub_avg:
	input:
		get_bigwig_rep_sub
	output:
		bigWigDir_sub_avg + "/{groupName}.subInput.avg.bw"
	message:
		"Making average bigWig files... [{wildcards.groupName}]"
	params:
		memory = "5G"
	shell:
		"""
		module load Cutlery/1.0
		makeBigWigAverage.sh -g {chrom_size} -m {params.memory} -o {output} {input}
		"""

def get_bigwig_rep_scaled_sub(wildcards):
	repL = samples.Name[samples.Group == wildcards.groupName].tolist()
	return map(lambda x: bigWigScaledDir_sub + "/" + x + ".subInput.bw", repL)

rule make_bigwig_scaled_sub_avg:
	input:
		get_bigwig_rep_scaled_sub
	output:
		bigWigScaledDir_sub_avg + "/{groupName}.scaled.subInput.avg.bw"
	message:
		"Making average bigWig files... [{wildcards.groupName}]"
	params:
		memory = "5G"
	shell:
		"""
		module load Cutlery/1.0
		makeBigWigAverage.sh -g {chrom_size} -m {params.memory} -o {output} {input}
		"""

def get_bigwig_rep_scaled_div(wildcards):
	repL = samples.Name[samples.Group == wildcards.groupName].tolist()
	return map(lambda x: bigWigScaledDir_div + "/" + x + ".divInput.bw", repL)

rule make_bigwig_scaled_div_avg:
	input:
		get_bigwig_rep_scaled_div
	output:
		bigWigScaledDir_div_avg + "/{groupName}.scaled.divInput.avg.bw"
	message:
		"Making average bigWig files... [{wildcards.groupName}]"
	params:
		memory = "5G"
	shell:
		"""
		module load Cutlery/1.0
		makeBigWigAverage.sh -g {chrom_size} -m {params.memory} -o {output} {input}
		"""
'''
