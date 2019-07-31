#!/usr/bin/env Rscript


#suppressPackageStartupMessages(library('geneplotter', quiet=TRUE))
suppressPackageStartupMessages(library('RColorBrewer', quiet=TRUE))
#suppressPackageStartupMessages(library('MASS', quiet=TRUE))
#suppressPackageStartupMessages(library('robustbase', quiet=TRUE))
suppressPackageStartupMessages(library('optparse', quiet=TRUE))
#suppressPackageStartupMessages(library('KernSmooth', quiet=TRUE))

source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))
# command line option handling
option_list <- list(
#	make_option(c("-g","--gtf"), default=NULL, help="gene annotation file"),
#	make_option(c("-r","--rRNA"), default=NULL, help="rRNA annotation file"),
#	make_option(c("-u","--update"), default=FALSE, action="store_true", help="Update existing record, default=FALSE"),
#	make_option(c("-s","--strand"), default=".", help="Strand requirement, + (same) / - (different) / . (both), default=."),
#	make_option(c("-f","--filtered"), default=FALSE, action="store_true", help="Strand requirement, + (same) / - (different) / . (both), default=.")

#	make_option(c("-y","--ylabel"), default="yLabel", help="Y-axis label"),
#	make_option(c("-t","--title"), default="Title", help="Main Title [default: Title]"),
#	make_option(c("-s","--size"), default="600,600", help="Comma-separated figure size, xSize,ySize"),
#	make_option(c("-o","--outfile"), default="barplot.png", help="Output file, with .png extension. [default: barplot.png]")
#	make_option(c("-h","--help"), default=FALSE, action="store_true", help="Print usage")
)
parser <- OptionParser(usage = "%prog [options] <list of alignment directory>", option_list=option_list,
			description = "For a give list of STAR alignment directories, print alignment stats")
arguments <- parse_args(parser, positional_arguments = TRUE)
if(length(arguments$args) == 0) {
	print_help(parser)
	stop("Requires alignment directory")
} else {
	srcDirL <- arguments$args
}

# Option handling
opt=arguments$options

#dataDirL=c("FGC0588_1_CGATGT", "FGC0588_1_GAGTGG", "FGC0588_1_GTGAAA")


if( FALSE ){
	gtf="../genes.uSymbol.gtf"
	srcDirL=c(
		"1.STAR/FGC1420_s_4_CGTACTAG_AGAGTAGA/",
		"1.STAR/FGC1420_s_4_CGTACTAG_GCGTAAGA/"
	)
	
	srcDirL=c(
		"ATAC_Day00_AAGAGGCA_AGAGTAGA_1",
		"ATAC_Day00_AAGAGGCA_AGAGTAGA_2",
		"ATAC_Day00_AAGAGGCA_CTCTCTAT_1",
		"ATAC_Day00_AAGAGGCA_CTCTCTAT_2",
		"ATAC_Day00_AAGAGGCA_GTAAGGAG_1"
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

for( srcDir in srcDirL ){
# srcDir=srcDirL[1]
#	write(sprintf("Checking STAR alignment: %s", srcDir), stderr())
	src=sprintf("%s/alignment.STARLog.final.out", srcDir)
	id=basename(srcDir)
#	stopifnot(file.exists(src))
	if( !file.exists(src) ){
		write(sprintf("%s\tNo information", id ), stderr())
		next
	}
#	stat = stat.template
	
	## STAR alignment STAT
	data = read.delim(sprintf("%s/alignment.STARLog.final.out", srcDir), header=FALSE, stringsAsFactor=FALSE)
	data = sapply(data, function(x) sub("[ |]+$", "", sub("^ *","", x)))
	data[,2] = sapply(data[,2], function(x) sub("%$","", x))	
	data = data[ grep(":$", data[,1], invert=TRUE), ]
	rownames(data) = data[,1]
#	data = data[unlist(fieldL),]
#	data = as.data.frame(data, stringsAsFactors=FALSE)
#	data[,2] = as.numeric(data[,2])


	stat = as.character(data[ unlist(fieldL), 2])
	write(paste(c(id,stat), collapse="\t"), stdout())
}




