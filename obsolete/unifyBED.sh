#!/usr/bin/env bash
# For a given set of BED files
# - sortBed
# - eliminate monoclonal tags leaving only one

tail -n +2 $1 | gawk '{OFS="\t"; print $1,$2,$3,"TAG",$5,$6}' | grep -v _ | sort -k1,1 -k2,2n -k6,6 | uniq
