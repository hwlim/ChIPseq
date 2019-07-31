###########################################################
# ChIP-seq pipeline using STAR aligner

###################
## For STAR
starIndex=/opt/STAR_Index/GENCODE.byLim/GRCh38_Gencode25/
outReadsUnmapped=None # None or Fastx
bamSort=TRUE	# if set, sorted bam is created


#######################
## Additional ChIPseq-specific options not allowing splicing
## if --aliginIntronMax is smaller than --aliginIntronMin, no denovo splicing is absolutely disabled.
## soft clipping
optStr_star="--runThreadN 4 --alignSJDBoverhangMin 999 --alignIntronMax 1 --alignMatesGapMax 2000 --outFilterMultimapNmax 1"

## bowtie1 style end-to-end mapping
optStr_star="--runThreadN 4 --alignSJDBoverhangMin 999 --alignIntronMax 1 --alignMatesGapMax 2000 --outFilterMultimapNmax 1 --alignEndsType EndToEnd --outFilterMismatchNoverLmax 0.05"


#########################
## General options
genome=hg38
## Directories
fqBase=../0.Fastq
alignDir=1.Align

## for chip.2.makeTagDir.sh
## If not NULL, chip.2.makeTagDir.sh will try to find alignment files in this directory
alignDir2=NULL


## Chromosome to select for downstream analysis. Comment it if N/A
chrRegex='chr[0-9XYM]*$'



