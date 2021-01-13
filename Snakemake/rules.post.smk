
rule count_spikein:
	input:
		fragDir + "/{sampleName}.frag.bed.gz"
	output:
		spikeinCntDir + "/{sampleName}.spikeCnt.txt"
	message:
		"Counting spikein tags... [{wildcards.sampleName}]"
	shell:
		"""
		module load Cutlery/1.0
		countSpikein.sh -p {spikePrefix} {input} > {output}
		"""

rule make_spikeintable:
	input:
		expand(spikeinCntDir + "/{sampleName}.spikeCnt.txt", sampleName=samples.Name.tolist())
	output:
		spikeinCntDir + "/spikein.txt"
	message:
		"Making spikein table..."
	shell:
		"""
		module load Cutlery/1.0
		makeSpikeCntTable.r -o {spikeinCntDir}/spikein {input}
		"""

## Draw a plot of fragment length distribution
## Automatically consider targer chromosomes only starting with "chr" excluding others such as "DM-chr"
rule get_fragLenHist:
	input:
		#dedupDir + "/{sampleName}.dedup.bam"
		fragDir + "/{sampleName}.frag.bed.gz"
	output:
		fragLenDir + "/{sampleName}.txt",
		fragLenDir + "/{sampleName}.png"
	message:
		"Checking fragment length... [{wildcards.sampleName}]"
	shell:
		"""
		module load Cutlery/1.0
		ngs.fragLenHist.r -o {fragLenDir}/{wildcards.sampleName} -n {wildcards.sampleName} {input}
		"""

## bigwig file: resized fragment in RPM scale
## Draw a plot of fragment length distribution
## Automatically consider targer chromosomes only starting with "chr" excluding others such as "DM-chr"
## and the same in other bigwig file rules
rule make_bigwig_ctr_rpm:
	input:
		fragDir_ctr + "/{sampleName}.frag.bed.gz"
	output:
		bigWigDir_ctr_RPM + "/{sampleName}.ctr.rpm.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	params:
		memory = "%dG" % (  cluster["make_bigwig"]["memory"]/1000 - 2 )
	shell:
		"""
		module load Cutlery/1.0
		cnr.fragToBigWig.sh -g {chrom_size} -m {params.memory} -o {output} {input}
		"""

## bigwig file: original fragment in RPM scale
rule make_bigwig_frag_rpm:
	input:
		fragDir + "/{sampleName}.frag.bed.gz"
	output:
		bigWigDir_frag_RPM + "/{sampleName}.frag.rpm.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	params:
		memory = "%dG" % (  cluster["make_bigwig"]["memory"]/1000 - 2 )
	shell:
		"""
		module load Cutlery/1.0
		cnr.fragToBigWig.sh -g {chrom_size} -m {params.memory} -o {output} {input}
		"""


## bigwig file: resized fragment in RPSM scale
rule make_bigwig_ctr_rpsm:
	input:
		bed = fragDir_ctr + "/{sampleName}.frag.bed.gz",
		spikeinCnt = spikeinCntDir + "/{sampleName}.spikeCnt.txt"
		#spikeinCnt = spikeinCntDir + "/spikein.txt"
	output:
		bigWigDir_ctr_RPSM + "/{sampleName}.ctr.rpsm.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	params:
		memory = "%dG" % ( cluster["make_bigwig"]["memory"]/1000 - 2 )
	shell:
		"""
		module load Cutlery/1.0
		scaleFactor=`cat {input.spikeinCnt} | gawk '{{ if($1=="ScaleFactor") print $2 }}'`
		if [ $scaleFactor == "" ];then
			echo -e "Error: empty scale factor" >&2
			exit 1
		fi
		cnr.fragToBigWig.sh -g {chrom_size} -m 5G -s $scaleFactor -o {output} {input.bed}
		"""

## bigwig file: original fragment in RPSM scale
rule make_bigwig_frag_rpsm:
	input:
		bed = fragDir + "/{sampleName}.frag.bed.gz",
		spikeinCnt = spikeinCntDir + "/{sampleName}.spikeCnt.txt"
		#spikeinCnt = spikeinCntDir + "/spikein.txt"
	output:
		bigWigDir_frag_RPSM + "/{sampleName}.frag.rpsm.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	params:
		memory = "%dG" % (  cluster["make_bigwig"]["memory"]/1000 - 2 )
	shell:
		"""
		module load Cutlery/1.0
		scaleFactor=`cat {input.spikeinCnt} | gawk '{{ if($1=="ScaleFactor") print $2 }}'`
		if [ $scaleFactor == "" ];then
			echo -e "Error: empty scale factor" >&2
			exit 1
		fi
		cnr.fragToBigWig.sh -g {chrom_size} -m 5G -s $scaleFactor -o {output} {input.bed}
		"""
#		scaleFactor=`cat {input.spikeinCnt} | gawk '$1=="'{wildcards.sampleName}'"' | gawk '{{ printf "%f", 100000/$3 }}'`



'''
rule make_bigwig_allfrag:
	input:
		fragDir + "/{sampleName}.frag.bed.gz"
	output:
		bigWigDirAllFrag + "/{sampleName}.allFrag.bw"
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
#	params:
#		memory = "5G"
	shell:
		"""
		module load Cutlery/1.0
		cnr.fragToBigWig.sh -g {chrom_size} -m 5G -o {output} {input}
		"""

rule make_bigwig_allfrag_rpsm:
	input:
		bed = fragDir + "/{sampleName}.frag.bed.gz",
		spikeinCnt = spikeinCntDir + "/spikein.txt"
	output:
		bigWigAllFrag_RPSM + "/{sampleName}.allFrag.rpsm.bw"
	message:
		"Making spike-in scaled allFrag bigWig files... [{wildcards.sampleName}]"
#	params:
#		memory = "5G"
	shell:
		"""
		module load Cutlery/1.0
		scaleFactor=`cat {input.spikeinCnt} | gawk '$1=="'{wildcards.sampleName}'"' | gawk '{{ printf "%f", 100000/$3 }}'`
		if [ $scaleFactor == "" ];then
			echo -e "Error: empty scale factor" >&2
			exit 1
		fi
		cnr.fragToBigWig.sh -g {chrom_size} -m 5G -s $scaleFactor -o {output} {input.bed}
		"""
'''


def get_input_name(sampleName):
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
		chip = bigWigDir_ctr_RPM + "/{sampleName}.ctr.rpm.bw",
		ctrl = lambda wildcards: bigWigDir_ctr_RPM + "/" + get_input_name(wildcards.sampleName) + ".ctr.rpm.bw"
	output:
		bigWigDir_ctr_RPM_sub + "/{sampleName}.ctr.rpm.subInput.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	#params:
	#	memory = "%dG" % (  cluster["make_bigwig_subtract"]["memory"]/1000 - 1 )
	shell:
		"""
		module load Cutlery/1.0
		bigWigSubtract.sh -g {chrom_size} -m 5G -t -1000 {output} {input.chip} {input.ctrl}
		"""

rule make_bigwig_frag_rpm_subinput:
	input:
		chip = bigWigDir_frag_RPM + "/{sampleName}.frag.rpm.bw",
		ctrl = lambda wildcards: bigWigDir_frag_RPM + "/" + get_input_name(wildcards.sampleName) + ".frag.rpm.bw"
	output:
		bigWigDir_frag_RPM_sub + "/{sampleName}.frag.rpm.subInput.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	#params:
	#	memory = "%dG" % (  cluster["make_bigwig_subtract"]["memory"]/1000 - 1 )
	shell:
		"""
		module load Cutlery/1.0
		bigWigSubtract.sh -g {chrom_size} -m 5G -t -1000 {output} {input.chip} {input.ctrl}
		"""


rule make_bigwig_ctr_rpsm_subinput:
	input:
		chip = bigWigDir_ctr_RPSM + "/{sampleName}.ctr.rpsm.bw",
		ctrl = lambda wildcards: bigWigDir_ctr_RPSM + "/" + get_input_name(wildcards.sampleName) + ".ctr.rpsm.bw"
	output:
		bigWigDir_ctr_RPSM_sub + "/{sampleName}.ctr.rpsm.subInput.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	#params:
	#	memory = "%dG" % (  cluster["make_bigwig_subtract"]["memory"]/1000 - 1 )
	shell:
		"""
		module load Cutlery/1.0
		bigWigSubtract.sh -g {chrom_size} -m 5G -t -1000 {output} {input.chip} {input.ctrl}
		"""

rule make_bigwig_frag_rpsm_subinput:
	input:
		chip = bigWigDir_frag_RPSM + "/{sampleName}.frag.rpsm.bw",
		ctrl = lambda wildcards: bigWigDir_frag_RPSM + "/" + get_input_name(wildcards.sampleName) + ".frag.rpsm.bw"
	output:
		bigWigDir_frag_RPSM_sub + "/{sampleName}.frag.rpsm.subInput.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	#params:
	#	memory = "%dG" % (  cluster["make_bigwig_subtract"]["memory"]/1000 - 1 )
	shell:
		"""
		module load Cutlery/1.0
		bigWigSubtract.sh -g {chrom_size} -m 5G -t -1000 {output} {input.chip} {input.ctrl}
		"""


'''
rule make_bigwig_allfrag_rpsm_subtract:
	input:
		chip = bigWigAllFrag_RPSM + "/{sampleName}.allFrag.rpsm.bw",
		ctrl = lambda wildcards: bigWigAllFrag_RPSM + "/" + get_input_name(wildcards.sampleName) + ".allFrag.rpsm.bw",
	output:
		bigWigAllFrag_RPSM_subInput + "/{sampleName}.allFrag.rpsm.subInput.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	#params:
	#	memory = "%dG" % (  cluster["make_bigwig_subtract"]["memory"]/1000 - 1 )
	shell:
		"""
		module load Cutlery/1.0
		bigWigSubtract.sh -g {chrom_size} -m 5G -t -1000 {output} {input.chip} {input.ctrl}
		"""


def get_bigwig_input(wildcards):
	# return ordered [ctrl , target] list.
	ctrlName = samples.Ctrl[samples.Name == wildcards.sampleName]
	ctrlName = ctrlName.tolist()[0]
	return [ bigWigDir + "/" + wildcards.sampleName + ".bw", bigWigDir + "/" + ctrlName + ".bw" ]

## RPM-scaled bigWig input-subtracted
rule make_bigwig_subtract:
	input:
		get_bigwig_input
	output:
		bigWigDir_sub + "/{sampleName}.subInput.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	#params:
	#	memory = "%dG" % (  cluster["make_bigwig_subtract"]["memory"]/1000 - 1 )
	shell:
		"""
		module load Cutlery/1.0
		bigWigSubtract.sh -g {chrom_size} -m 5G -t -1000 {output} {input}
		"""


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
'''

def get_hetchr_ctrl(sampleName):
	#print(sampleName)
	ctrlName = samples.Ctrl[samples.Name == sampleName]
	ctrlName = ctrlName.tolist()[0]
	return ctrlName
#	return [ fragDir + "/" + sampleName + ".frag.bed.gz", fragDir + "/" + ctrlName + ".frag.bed.gz" ]


rule call_peak_hetchr:
	input:
		chip = fragDir + "/{sampleName}.frag.bed.gz",
		ctrl = lambda wildcards: fragDir + "/" + get_hetchr_ctrl(wildcards.sampleName) + ".frag.bed.gz",
		spikeChip = spikeinCntDir + "/{sampleName}.spikeCnt.txt",
		spikeCtrl = lambda wildcards: spikeinCntDir + "/" + get_hetchr_ctrl(wildcards.sampleName) + ".spikeCnt.txt"
	output:
		peakDir + "/{sampleName}." + peakSuffix + ".bed",
		peakDir + "/{sampleName}." + peakSuffix + ".txt.gz"
	params:
		outPrefix = peakDir + "/{sampleName}." + peakSuffix,

	message:
		"Peak calling for heterochromatin... [{wildcards.sampleName}]"
	shell:
		"""
		module load ChIPseq
		spikeChip=`cat {input.spikeChip} | grep ^Spikein | cut -f 2`
		spikeCtrl=`cat {input.spikeCtrl} | grep ^Spikein | cut -f 2`
		if [ "$spikeChip" == "" ] || [ "$spikeCtrl" == "" ];then
			echo -e "Error: empty spikein factor" >&2
			exit 1
		fi
		spikein=`echo -e "${{spikeChip}}\t${{spikeCtrl}}" | gawk '{{ printf "%f", $2 / $1 }}'`
		hetChr.call.sh -o {params.outPrefix} -g {chrom_size} -w {peakWindow} -s {peakStep} -a {peakAlpha} -k $spikein -t {peakFC} {input.chip} {input.ctrl}
		"""

rule call_peak_hetchr2:
	input:
		chip = fragDir + "/{sampleName}.frag.bed.gz",
		ctrl = lambda wildcards: fragDir + "/" + get_hetchr_ctrl(wildcards.sampleName) + ".frag.bed.gz",
		spikeChip = spikeinCntDir + "/{sampleName}.spikeCnt.txt",
		spikeCtrl = lambda wildcards: spikeinCntDir + "/" + get_hetchr_ctrl(wildcards.sampleName) + ".spikeCnt.txt"
	output:
		peakDir_homer + "/{sampleName}." + peakSuffix_homer + ".exBL.bed",
	params:
		outPrefix = peakDir_homer + "/{sampleName}." + peakSuffix_homer,

	message:
		"Peak calling for heterochromatin... [{wildcards.sampleName}]"
	shell:
		"""
		module load ChIPseq
		spikeChip=`cat {input.spikeChip} | grep ^Spikein | cut -f 2`
		spikeCtrl=`cat {input.spikeCtrl} | grep ^Spikein | cut -f 2`
		if [ "$spikeChip" == "" ] || [ "$spikeCtrl" == "" ];then
			echo -e "Error: empty spikein factor" >&2
			exit 1
		fi
		spikein=`echo -e "${{spikeChip}}\t${{spikeCtrl}}" | gawk '{{ printf "%f", $1 / $2 }}'`
		hetChr.call.homer.sh -o {params.outPrefix} -s {peakWindow_homer} -d {minDist_homer} -k $spikein -f {peakFC_homer} {input.chip} {input.ctrl}
		"""


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
		if [ $scaleFactor == "" ];then
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