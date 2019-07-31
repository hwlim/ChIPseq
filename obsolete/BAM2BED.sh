#!/usr/bin/env bash

#python ../scripts/S0.BAM2BED.Cluster.py $1 $2

# not yet tested below
#echo $1, $2
bamToBed -i $1 | gawk '{OFS="\t"; if(strtonum($5)>'$3') print $1,$2,$3,$4,$5,$6}' > $2 
