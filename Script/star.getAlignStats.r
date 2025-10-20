#!/usr/bin/env Rscript


#suppressPackageStartupMessages(library('geneplotter', quiet=TRUE))
suppressPackageStartupMessages(library('RColorBrewer', quiet=TRUE))
#suppressPackageStartupMessages(library('MASS', quiet=TRUE))
#suppressPackageStartupMessages(library('robustbase', quiet=TRUE))
suppressPackageStartupMessages(library('optparse', quiet=TRUE))
#suppressPackageStartupMessages(library('KernSmooth', quiet=TRUE))
library(ggplot2)
library(cowplot)
library(reshape2)

source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))
# command line option handling
option_list <- list(
	make_option(c("-o","--outPrefix"), default=NULL, help="Output prefix for alignment statistics & plot files, default=NULL
		Plot size is automatically decided by the sample number")
#	make_option(c("-r","--rRNA"), default=NULL, help="rRNA annotation file"),
#	make_option(c("-u","--update"), default=FALSE, action="store_true", help="Update existing record, default=FALSE"),
#	make_option(c("-s","--strand"), default=".", help="Strand requirement, + (same) / - (different) / . (both), default=."),
#	make_option(c("-f","--filtered"), default=FALSE, action="store_true", help="Strand requirement, + (same) / - (different) / . (both), default=.")

#	make_option(c("-y","--ylabel"), default="yLabel", help="Y-axis label"),
#	make_option(c("-t","--title"), default="Title", help="Main Title [default: Title]"),
#	make_option(c("-s","--size"), default="12,600", help="Comma-separated figure size, xSize,ySize"),
#	make_option(c("-o","--outfile"), default="barplot.png", help="Output file, with .png extension. [default: barplot.png]")
#	make_option(c("-h","--help"), default=FALSE, action="store_true", help="Print usage")
)
parser <- OptionParser(usage = "%prog [options] <list of alignment directory>", option_list=option_list,
			description = "For a give list of STAR alignment directories, print alignment stats to stdio
If outPrefix is specified, <outPrefix>.txt and <outPrefix>.png files are also created.")
arguments <- parse_args(parser, positional_arguments = TRUE)
if(length(arguments$args) == 0) {
	print_help(parser)
	stop("Requires alignment directory")
} else {
	srcDirL <- arguments$args
}

# Option handling
opt=arguments$options
outPrefix = opt$outPrefix
#dataDirL=c("FGC0588_1_CGATGT", "FGC0588_1_GAGTGG", "FGC0588_1_GTGAAA")


if( FALSE ){
	srcDirL=c(
		"../1.1.Align/DE_H3K9me3_Control_A6",
		"../1.1.Align/DE_H3K9me3_Control_H5",
		"../1.1.Align/DE_H3K9me3_FOXA123KD_A6"
	)

}



fieldL=list()
fieldL[["TTC"]]="Number of input reads"
fieldL[["AlnUniq"]]="Uniquely mapped reads number"
fieldL[["AlnMulti"]]="Number of reads mapped to multiple loci"

fieldL[["AlnUniq%"]]="Uniquely mapped reads %"
fieldL[["AlnMulti%"]]="% of reads mapped to multiple loci"
fieldL[["FailMulti%"]]="% of reads mapped to too many loci"
fieldL[["FailMisMatch%"]]="% of reads unmapped: too many mismatches"
fieldL[["FailTooShort%"]]="% of reads unmapped: too short"
fieldL[["FailOther%"]]="% of reads unmapped: other"


write(paste(c("Sample", names(fieldL)), collapse="\t"), stdout())

df.stat = NULL
for( srcDir in srcDirL ){
# srcDir=srcDirL[1]
#	write(sprintf("Checking STAR alignment: %s", srcDir), stderr())
	# srcDir="hPSC_Input_Control_A6"
	src=sprintf("%s/%s", srcDir, list.files(srcDir, pattern=".final.out")[1])
	#src=sprintf("%s/alignment.STARLog.final.out", srcDir)
	id=basename(srcDir)
#	stopifnot(file.exists(src))
	if( !file.exists(src) ){
		write(sprintf("%s\tNo information", id ), stderr())
		next
	}
#	stat = stat.template
	
	## STAR alignment STAT
	data = read.delim(src, header=FALSE, stringsAsFactor=FALSE)
	data = sapply(data, function(x) sub("[ |]+$", "", sub("^ *","", x)))
	data[,2] = sapply(data[,2], function(x) sub("%$","", x))	
	data = data[ grep(":$", data[,1], invert=TRUE), ]
	rownames(data) = data[,1]
#	data = data[unlist(fieldL),]
#	data = as.data.frame(data, stringsAsFactors=FALSE)
#	data[,2] = as.numeric(data[,2])


	stat = as.character(data[ unlist(fieldL), 2])
	write(paste(c(id,stat), collapse="\t"), stdout())

	if(is.null(df.stat)){
		df.stat = as.numeric(stat)
	}else{
		df.stat = rbind(df.stat, as.numeric(stat))
	}
}
df.stat = data.frame(df.stat)
rownames(df.stat) = basename(srcDirL)
colnames(df.stat) = names(fieldL)
df.frac = df.stat[,grep("%$", names(df.stat))]
colnames(df.frac) = sub("%$","", colnames(df.frac))


if( !is.null(outPrefix) ){
	## Single barplot
	df.melt = melt(data.frame(Sample=rownames(df.frac), df.frac))
	colnames(df.melt) = c("Sample","Alignment","Percent")
	df.melt$Alignment = factor(df.melt$Alignment, levels=colnames(df.frac))
	g = ggplot(df.melt, aes(fill=Alignment, y=Percent, x=Sample)) +
		geom_bar(position="stack", stat="identity") +
		coord_flip()
	ggsave(sprintf("%s.single.pdf", outPrefix),  g, width=9, height=max(4, 2 + 0.2*length(srcDirL)))


	## Multipple barplots stratified by alignment type / failure
	df.stat = data.frame(Sample = rownames(df.stat), df.stat)
	df.stat$Sample = factor(df.stat$Sample, levels=df.stat$Sample)

	g.ttc = ggplot(data=df.stat, aes(x=Sample, y=TTC)) +
		geom_bar(stat="identity", fill="steelblue")+
		theme_minimal() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
		labs(title="Total Reads", y="Count")

	g.uniq = ggplot(data=df.stat, aes(x=Sample, y=AlnUniq)) +
		geom_bar(stat="identity", fill="steelblue")+
		theme_minimal() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
		labs(title="Uniquely Aligneds", y="Count")

	g.uniqRate = ggplot(data=df.stat, aes(x=Sample, y=AlnUniq.)) +
		geom_bar(stat="identity", fill="steelblue")+
		theme_minimal() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
		labs(title="Uniquely Aligneds (%)", y="%")

	g.failMulti = ggplot(data=df.stat, aes(x=Sample, y=FailMulti.)) +
		geom_bar(stat="identity", fill="steelblue")+
		theme_minimal() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
		labs(title="Fail by Multi Align (%)", y="%")

	g.failShort = ggplot(data=df.stat, aes(x=Sample, y=FailTooShort.)) +
		geom_bar(stat="identity", fill="steelblue")+
		theme_minimal() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
		labs(title="Fail by Too Short (%)", y="%")

	g.failOther = ggplot(data=df.stat, aes(x=Sample, y=FailOther.)) +
		geom_bar(stat="identity", fill="steelblue")+
		theme_minimal() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
		labs(title="Fail Other (%)", y="%")

	g = plot_grid( g.ttc, g.uniq, g.uniqRate, g.failMulti, g.failShort, g.failOther, ncol=3)

	if(length(srcDirL) < 10){
		w.unit = 0.2
	}else{
		w.unit = 0.1
	}
	ggsave(sprintf("%s.split.pdf", outPrefix), g, height=6, width=(1+w.unit*length(srcDirL))*3)
	write.table(df.stat, sprintf("%s.txt", outPrefix), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
}

