

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
		module load CnR/1.0
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

rule count_spikein:
	input:
		fragDir + "/{sampleName}.frag.bed.gz"
	output:
		spikeinCntDir + "/{sampleName}.spikeCnt.txt"
	message:
		"Counting spikein tags... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
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
		module load CnR/1.0
		makeSpikeCntTable.r -o {spikeinCntDir}/spikein {input}
		"""

rule get_fragLenHist:
	input:
		dedupDir + "/{sampleName}.dedup.bam"
	output:
		fragLenDir + "/{sampleName}.dist.txt",
		fragLenDir + "/{sampleName}.dist.png"
	message:
		"Checking fragment length... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		ngs.fragLenHist.r -o {fragLenDir}/{wildcards.sampleName} {input}
		"""

## RPM-scaled bigWig
rule make_bigwig:
	input:
		fragDir + "/{sampleName}.frag.bed.gz"
	output:
		bigWigDir + "/{sampleName}.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	params:
		memory = "%dG" % (  cluster["make_bigwig"]["memory"]/1000 - 1 )
	shell:
		"""
		module load CnR/1.0
		cnr.bedToBigWig.sh -g {chrom_size} -m {params.memory} -o {output} {input}
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
		module load CnR/1.0
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
		module load CnR/1.0
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
		module load CnR/1.0
		cnr.makeHomerDir.sh -c {chrRegexTarget} -o {output} -n {params.name} {input}
		"""

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
'''
rule call_peak_hetchr2:
	input:
		chip = fragDir + "/{sampleName}.frag.bed.gz",
		ctrl = lambda wildcards: fragDir + "/" + get_hetchr_ctrl(wildcards.sampleName) + ".frag.bed.gz",
		spikeChip = spikeinCntDir + "/{sampleName}.spikeCnt.txt",
		spikeCtrl = lambda wildcards: spikeinCntDir + "/" + get_hetchr_ctrl(wildcards.sampleName) + ".spikeCnt.txt"
	output:
		peakDir2 + "/{sampleName}." + peakSuffix + ".exBL.bed",
	params:
		outPrefix = peakDir2 + "/{sampleName}." + peakSuffix,

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
		hetChr.call.homer.sh -o {params.outPrefix} -s {peakWindow} -k $spikein -t {peakFC} {input.chip} {input.ctrl}
		"""
		#hetChr.call.sh -o {params.outPrefix} -g {chrom_size} -w {peakWindow} -s {peakStep} -a {peakAlpha} -k $spikein -t {peakFC} {input.chip} {input.ctrl}
...

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
		module load CnR/1.0
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
		module load CnR/1.0
		scaleFactor=`cat {input.spikeinCnt} | gawk '$1=="'{wildcards.sampleName}'"' | cut -f 6`
		if [ $scaleFactor == "" ];then
			echo -e "Error: empty scale factor" >&2
			exit 1
		fi
		cnr.bedToBigWig.sh -g {chrom_size} -m 5G -s $scaleFactor -o {output} {input.bed}
		"""
#		cnr.bedToBigWig.sh -g {chrom_size} -m {params.memory} -s {params.scaleFactor} -o {output} {input.bed}

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
		module load CnR/1.0
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
		module load CnR/1.0
		bigWigDivide.sh -g {chrom_size} -m 5G -s log -a 1 -o {output} {input}
		"""




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
		module load CnR/1.0
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
		module load CnR/1.0
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
		module load CnR/1.0
		makeBigWigAverage.sh -g {chrom_size} -m {params.memory} -o {output} {input}
		"""
