#!/usr/bin/env bash

## Args
chromSize=$1
peakMask=$2
targetBam=$3
inputBam=$4
outDir=$5
eddConfig=$6

## Create temp files
tmpLogRatioSorted=$(mktemp "${outDir}/EDD_log_ratio_sorted.XXXXXX.bedgraph")
tmpBinScoreSorted=$(mktemp "${outDir}/EDD_bin_score_sorted.XXXXXX.bedgraph")
tmpLogRatioFixed=$(mktemp "${outDir}/EDD_log_ratio_fixed.XXXXXX.bedgraph")
tmpBinScoreFixed=$(mktemp "${outDir}/EDD_bin_score_fixed.XXXXXX.bedgraph")
tmpLogRatioZeros=$(mktemp "${outDir}/EDD_log_ratio_zeros.XXXXXX.bedgraph")
tmpBinScoreZeros=$(mktemp "${outDir}/EDD_bin_score_zeros.XXXXXX.bedgraph")

## Set up cleanup trap
trap "rm -f $tmpLogRatioSorted $tmpBinScoreSorted $tmpLogRatioFixed $tmpBinScoreFixed $tmpLogRatioZeros $tmpBinScoreZeros" EXIT

## run EDD
edd --config-file $eddConfig --write-bin-scores --write-log-ratios $chromSize $peakMask $targetBam $inputBam $outDir

## sort output bedgraph files
sort -k1,1 -k2,2n ${outDir}/EDD_log_ratio_1KB.bedgraph > $tmpLogRatioSorted
sort -k1,1 -k2,2n ${outDir}/EDD_bin_score.bedgraph > $tmpBinScoreSorted

## remove coordinates that are out of bounds
sed 's/[[:space:]]*$//' $tmpLogRatioSorted \
    | gawk 'BEGIN{OFS="\t"} 
        NR==FNR {chrom_size[$1]=$2; next}
        {
        if ($3 > chrom_size[$1]) $3 = chrom_size[$1]
        if ($2 < chrom_size[$1]) print
        }' $chromSize - > $tmpLogRatioFixed

sed 's/[[:space:]]*$//' $tmpBinScoreSorted \
    | gawk 'BEGIN{OFS="\t"} 
        NR==FNR {chrom_size[$1]=$2; next}
        {
        if ($3 > chrom_size[$1]) $3 = chrom_size[$1]
        if ($2 < chrom_size[$1]) print
        }' $chromSize - > $tmpBinScoreFixed

## add zeros to 4th column; bedgraphtobigwig can't handle empty values
gawk 'BEGIN {FS=OFS="\t"} NF==3 {$4=0} 1' $tmpLogRatioFixed > $tmpLogRatioZeros
gawk 'BEGIN {FS=OFS="\t"} NF==3 {$4=0} 1' $tmpBinScoreFixed > $tmpBinScoreZeros

bedGraphToBigWig $tmpLogRatioZeros $chromSize ${outDir}/EDD_log_ratio_1KB.bw
bedGraphToBigWig $tmpBinScoreZeros $chromSize ${outDir}/EDD_bin_score.bw

mv ${outDir}/EDD_log_ratio_1KB.bw ${outDir}/igv.log_ratio.bw
mv ${outDir}/EDD_bin_score.bw ${outDir}/igv.bin_score.bw
mv ${outDir}/EDD_peaks.bed ${outDir}/peak.exBL.bed 

rm ${outDir}/EDD_log_ratio_1KB.bedgraph ${outDir}/EDD_bin_score.bedgraph
