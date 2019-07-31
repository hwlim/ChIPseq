#!/usr/bin/env bash

# Written by Hee Woong Lim, 11/29/2011
# FGC bed file --> 6 column sorted, change read name, optionally read length modification

source $COMMON_LIB_BASE/commonBash.sh
trap "kill_child_processes 1 $$; exit 1" INT 
#set -e
#set -u

function printUsage {
	echo "Usage: normalizeAlign.sh (options) [alignFile]" >&2
	echo -e "Normalize align file (bam or bed, or .gz) and print out to stdio" >&2
	echo -e "ASSUMES TAB-deliminated FORMAT!" >&2
	echo -e "\tchr start end name count direc" >&2
	echo -e "Options:" >&2
	echo -e "\t-l #: Change the read length by extending or shrink at 3' end" >&2
	#echo -e "\t-o <dir>: Output directory name" >&2
	#echo -e "\t-4: if set, assume 4 field bed format, [chr,start,end,direc]" >&2
	#echo -e "\t-1: if set, ignore previous tag count and just set 1" >&2
	#echo -e "\t-u: Leave only one tag at each position, replace read name with 'read'" >&2
	echo -e "\t-d: delimeter, default=<tab>" >&2
	echo -e "\t-f <#,#,#,#,#>: chr,start,end,direc,count. 'count' field can be omitted like 1,2,3,6," >&2
	echo -e "\t-M : exclude chrM" >&2 
	echo -e "\t-s : sortBed" >&2 
	#echo -e "\t-r <#>: Down sampling rate, [0.1, 0.2, 0.3, 0.4, 0.5]" >&2
	echo -e "\t-c : Make compact file by collecting same tag positions to tag count" >&2
	echo -e "\t-1 : Input records are 1-based, need to subtract 1 from start coordinate (default: zero-based like BED)\n\t\t<NOT IMPLEMENTED YET>" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi

#optStr=""
#outDir="."
rdLen=0
#isTC1=0
#isBed4=0
#toUniq=0
#tbp=0
delimeter="\t"
field=""
toCompact=0
toExChrM=0
toSort=0
oneBase=0
while getopts ":l:f:d:Msc1" opt; do
	case $opt in
		l)
			rdLen=$OPTARG
			;;
		d)
			delimeter="$OPTARG"
			;;
		f)
			field=$OPTARG
			;;
		M)
			toExChrM=1
			;;
		s)
			toSort=1
			;;
		c)
			toCompact=1
			;;
		1)
			oneBase=1
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			printUsage
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			printUsage
			exit 1
			;;
	esac
done


#echo "$optStr"
shift $((OPTIND-1))

assertFileExist $@
alignL=($@)
align=${alignL[0]}

nF=`echo $field | grep -o "," | wc -l`
if [ $nF -ne 4 ];then
	echo "Invalid field format: %field" >&2
	printUsage
	exit 1
fi

catCmd="cat"
if [ "${align##*.}" == "gz" ];then
	catCmd="zcat"
elif [ "${align##*.}" == "bam" ];then
	catCmd="bamToBed -i"
#	catCmd="bam2bed.sh"
fi

if [ $toExChrM -eq 1 ];then
	regStr="^chr[0-9]*[XY]*[[:space:]]"
else
	regStr="^chr[0-9]*[XYM]*[[:space:]]"
fi

lastCmd="grep ${regStr}"
if [ $toSort -eq 1 ];then
       lastCmd="${lastCmd} | sort -S 5G -k1,1 -k2,2n -k3,3n -k6,6"
fi
if [ $toCompact -eq 1 ];then
	lastCmd="${lastCmd} | rd2compact.sh"
fi

echo "Processing ${alignL[@]}" >&2
echo "chr/start/end/direc/cnt: $field" >&2
echo "$lastCmd" >&2

function normalize {
	${catCmd} ${alignL[@]} \
	| gawk -F "$delimeter" 'BEGIN{
			fStr="'${field}'";
			rdLen='${rdLen}';
			split(fStr, fL, ",");
			chrInd=fL[1]; startInd=fL[2]; endInd=fL[3]; direcInd=fL[4]; cntInd=fL[5];
		}
		{
			chr=$(chrInd); start=$(startInd); end=$(endInd); direc=$(direcInd);
			if(cntInd==""){ cnt=1 }else{ cnt=$(cntInd) }
			if(rdLen>0){
				if(direc=="+"){
					end=start+rdLen;
				}else{
					start=end-rdLen;
				}
				if(start<=0){next}
			}
			printf "%s\t%d\t%d\t.\t%d\t%s\n", chr, start, end, cnt, direc;
		}'
}

lastCmd="grep ${regStr}"
if [ $toSort -eq 1 ] && [ $toCompact -eq 1 ];then
	normalize \
		| grep ${regStr}\
		| sort -S 5G -k1,1 -k2,2n -k3,3n -k6,6 \
		| rd2compact.sh
elif [ $toSort -eq 1 ] && [ $toCompact -eq 0 ];then
	normalize \
		| grep ${regStr}\
		| sort -S 5G -k1,1 -k2,2n -k3,3n -k6,6 
elif [ $toSort -eq 0 ] && [ $toCompact -eq 1 ];then
	normalize \
		| grep ${regStr}\
		| rd2compact.sh
else 
	normalize \
		| grep ${regStr}
fi

