#!/usr/bin/env Rscript

# R script to draw catter plot
# Usage:
#	drawScatter.r (option) [dataFile]
#	-x <xLabel>
#	-y <yLabel>
#	-t <mainTitle>
#	-o <outFile>
#	-l <-xlim,+xlim,-ylim,+ylim>, comma separated axis limit
#	-f : comma-separated field numbers for x/y axis and category
# Data file format:
#	

#suppressPackageStartupMessages(library('geneplotter', quiet=TRUE))
suppressPackageStartupMessages(library('RColorBrewer', quiet=TRUE))
#suppressPackageStartupMessages(library('MASS', quiet=TRUE))
suppressPackageStartupMessages(library('ggplot2', quiet=TRUE))
suppressPackageStartupMessages(library('optparse', quiet=TRUE))
#suppressPackageStartupMessages(library('KernSmooth', quiet=TRUE))
source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))

# command line option handling
option_list <- list(
	make_option(c("-t","--title"), default=NULL, help="Main Title [default: source file name]"),
	make_option(c("-o","--outPrefix"), default="hist.logfc", help="Output file, with .png extension. [default: hist.logfc]")
)
parser <- OptionParser(usage = "%prog [options] dataFile", option_list=option_list,
	description = "Description:
	For a give list of one-column text files, perform Venn diagram analysis.
Output:
	- <outPrefix>.<pdf,png?: Density plot of log2FC" )
arguments <- parse_args(parser, positional_arguments = TRUE)
if(length(arguments$args) == 0) {
	print_help(parser)
	stop("Error: Requires a data file")
} else {
	src <- arguments$args[1]
}

# Option handling
opt=arguments$options

title=opt$title
outPrefix=opt$outPrefix

assertFileExist(src)

if(is.null(title)) title=src
desDir=dirname(outPrefix)
system(sprintf("mkdir -p %s", desDir))

data = read.delim(src, header=FALSE)
colnames(data) = c("Chr","Start","End","log2FC")

g = ggplot(data, aes(x=log2FC))+
	geom_density() +
	ggtitle(title) +
	xlab("log2FC") + ylab("Density") +
	theme_bw()
ggsave(sprintf("%s.pdf", outPrefix), g, width=5, heigh=3)
ggsave(sprintf("%s.png", outPrefix), g, width=5, heigh=3)