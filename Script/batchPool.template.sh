#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh

configL=(
#Liver_CEBPb_GRdim_6am_pool/TSV		"mlds_500_5_CGATGT/TSVu mlds_500_5_TGACCA/TSVu"
#Liver_CEBPb_GRdim_6pm_pool/TSV		"mlds_500_5_ACAGTG/TSVu mlds_500_5_GCCAAT/TSVu"
#Liver_CEBPb_GRdim_pred_6am_pool/TSV	"mlds_500_5_CAGATC/TSVu mlds_500_5_CTTGTA/TSVu"
#Liver_CEBPb_WT_6am_pool/TSV		"mlds_500_4_CGATGT/TSVu mlds_500_4_TGACCA/TSVu"
#Liver_CEBPb_WT_6pm_pool/TSV		"mlds_500_4_ACAGTG/TSVu mlds_500_4_GCCAAT/TSVu"
#Liver_CEBPb_WT_pred_6am_pool/TSV	"mlds_500_4_CAGATC/TSVu mlds_500_4_CTTGTA/TSVu"

#Liver_GR_GRdim_6am_pool/TSV		"mlds_500_7_CGATGT/TSVu mlds_500_7_TGACCA/TSVu"
#Liver_GR_GRdim_6pm_pool/TSV		"mlds_500_7_ACAGTG/TSVu mlds_500_7_GCCAAT/TSVu"
#Liver_GR_GRdim_pred_6am_pool/TSV	"mlds_500_7_CAGATC/TSVu mlds_500_7_CTTGTA/TSVu"
#Liver_GR_WT_6am_pool/TSV		"mlds_500_6_CGATGT/TSVu mlds_605_6_ACAGTG/TSVu" 	#mlds_500_6_TGACCA/TSVu
#Liver_GR_WT_6pm_pool/TSV		"mlds_500_6_ACAGTG/TSVu mlds_500_6_GCCAAT/TSVu"
#Liver_GR_WT_pred_6am_pool/TSV		"mlds_500_6_CAGATC/TSVu mlds_500_6_CTTGTA/TSVu"

#Liver_GR_GRdim_6am_pool/TSVrd.15M	"mlds_500_7_CGATGT/TSVrd.15M mlds_500_7_TGACCA/TSVrd.15M"
#Liver_GR_GRdim_6pm_pool/TSVrd.15M	"mlds_500_7_ACAGTG/TSVrd.15M mlds_500_7_GCCAAT/TSVrd.15M"
#Liver_GR_GRdim_pred_6am_pool/TSVrd.15M	"mlds_500_7_CAGATC/TSVrd.15M mlds_500_7_CTTGTA/TSVrd.15M"
Liver_GR_WT_6am_pool/TSVrd.15M		"mlds_500_6_CGATGT/TSVrd.15M mlds_605_6_ACAGTG/TSVrd.15M"
#Liver_GR_WT_6pm_pool/TSVrd.15M		"mlds_500_6_ACAGTG/TSVrd.15M mlds_500_6_GCCAAT/TSVrd.15M"
#Liver_GR_WT_pred_6am_pool/TSVrd.15M	"mlds_500_6_CAGATC/TSVrd.15M mlds_500_6_CTTGTA/TSVrd.15M"

#Liver_CEBPb_GRdim_6am_pool/TSVrd.15M		"mlds_500_5_CGATGT/TSVrd.15M mlds_500_5_TGACCA/TSVrd.15M"
#Liver_CEBPb_GRdim_6pm_pool/TSVrd.15M		"mlds_500_5_ACAGTG/TSVrd.15M mlds_500_5_GCCAAT/TSVrd.15M"
#Liver_CEBPb_GRdim_pred_6am_pool/TSVrd.15M	"mlds_500_5_CAGATC/TSVrd.15M mlds_500_5_CTTGTA/TSVrd.15M"
#Liver_CEBPb_WT_6am_pool/TSVrd.15M		"mlds_500_4_CGATGT/TSVrd.15M mlds_500_4_TGACCA/TSVrd.15M"
#Liver_CEBPb_WT_6pm_pool/TSVrd.15M		"mlds_500_4_ACAGTG/TSVrd.15M mlds_500_4_GCCAAT/TSVrd.15M"
#Liver_CEBPb_WT_pred_6am_pool/TSVrd.15M		"mlds_500_4_CAGATC/TSVrd.15M mlds_500_4_CTTGTA/TSVrd.15M"

#Liver_GR_GRdim_6am_pool/TSV.15M		"mlds_500_7_CGATGT/TSVu.15M mlds_500_7_TGACCA/TSVu.15M"
#Liver_GR_GRdim_6pm_pool/TSV.15M		"mlds_500_7_ACAGTG/TSVu.15M mlds_500_7_GCCAAT/TSVu.15M"
#Liver_GR_GRdim_pred_6am_pool/TSV.15M	"mlds_500_7_CAGATC/TSVu.15M mlds_500_7_CTTGTA/TSVu.15M"
Liver_GR_WT_6am_pool/TSV.15M		"mlds_500_6_CGATGT/TSVu.15M mlds_605_6_ACAGTG/TSVu.15M"
#Liver_GR_WT_6pm_pool/TSV.15M		"mlds_500_6_ACAGTG/TSVu.15M mlds_500_6_GCCAAT/TSVu.15M"
#Liver_GR_WT_pred_6am_pool/TSV.15M	"mlds_500_6_CAGATC/TSVu.15M mlds_500_6_CTTGTA/TSVu.15M"

#Liver_CEBPb_GRdim_6am_pool/TSV.15M		"mlds_500_5_CGATGT/TSVu.15M mlds_500_5_TGACCA/TSVu.15M"
#Liver_CEBPb_GRdim_6pm_pool/TSV.15M		"mlds_500_5_ACAGTG/TSVu.15M mlds_500_5_GCCAAT/TSVu.15M"
#Liver_CEBPb_GRdim_pred_6am_pool/TSV.15M		"mlds_500_5_CAGATC/TSVu.15M mlds_500_5_CTTGTA/TSVu.15M"
#Liver_CEBPb_WT_6am_pool/TSV.15M			"mlds_500_4_CGATGT/TSVu.15M mlds_500_4_TGACCA/TSVu.15M"
#Liver_CEBPb_WT_6pm_pool/TSV.15M			"mlds_500_4_ACAGTG/TSVu.15M mlds_500_4_GCCAAT/TSVu.15M"
#Liver_CEBPb_WT_pred_6am_pool/TSV.15M		"mlds_500_4_CAGATC/TSVu.15M mlds_500_4_CTTGTA/TSVu.15M"
)

configL=(
Liver_Pol2_GRdim_6am_pool/TSVrd.15M	"mlds_629_8_CGATGT/TSVrd.15M mlds_629_8_TGACCA/TSVrd.15M"
Liver_Pol2_GRdim_6pm_pool/TSVrd.15M	"mlds_629_8_ACAGTG/TSVrd.15M mlds_629_8_GCCAAT/TSVrd.15M"
Liver_Pol2_GRdim_pred_6am_pool/TSVrd.15M	"mlds_629_8_CAGATC/TSVrd.15M mlds_629_8_CTTGTA/TSVrd.15M"
Liver_Pol2_WT_6am_pool/TSVrd.15M	"mlds_646_3_CGATGT/TSVrd.15M mlds_646_3_TGACCA/TSVrd.15M"
Liver_Pol2_WT_6pm_pool/TSVrd.15M	"mlds_646_3_ACAGTG/TSVrd.15M mlds_646_3_GCCAAT/TSVrd.15M"
Liver_Pol2_WT_pred_6am_pool/TSVrd.15M	"mlds_646_3_CAGATC/TSVrd.15M mlds_646_3_CTTGTA/TSVrd.15M"

Liver_Pol2_GRdim_6am_pool/TSV.15M	"mlds_629_8_CGATGT/TSVu.15M mlds_629_8_TGACCA/TSVu.15M"
Liver_Pol2_GRdim_6pm_pool/TSV.15M	"mlds_629_8_ACAGTG/TSVu.15M mlds_629_8_GCCAAT/TSVu.15M"
Liver_Pol2_GRdim_pred_6am_pool/TSV.15M	"mlds_629_8_CAGATC/TSVu.15M mlds_629_8_CTTGTA/TSVu.15M"
Liver_Pol2_WT_6am_pool/TSV.15M		"mlds_646_3_CGATGT/TSVu.15M mlds_646_3_TGACCA/TSVu.15M"
Liver_Pol2_WT_6pm_pool/TSV.15M		"mlds_646_3_ACAGTG/TSVu.15M mlds_646_3_GCCAAT/TSVu.15M"
Liver_Pol2_WT_pred_6am_pool/TSV.15M	"mlds_646_3_CAGATC/TSVu.15M mlds_646_3_CTTGTA/TSVu.15M"
)


configL=(
Liver_CEBPb_WT_6am_exo_pool/TSVrd	"mlds_646_4_CGATGT/TSVrd mlds_646_4_TGACCA/TSVrd"
Liver_CEBPb_WT_6pm_exo_pool/TSVrd	"mlds_646_4_CTTGTA/TSVrd mlds_646_4_GCCAAT/TSVrd"
Liver_GR_GRdim_6pm_exo_pool/TSVrd	"mlds_605_6_GCCAAT/TSVrd mlds_605_6_CTTGTA/TSVrd"
Liver_GR_WT_6am_exo_pool/TSVrd		"mlds_633_7_CGATGT/TSVrd mlds_633_7_TGACCA/TSVrd"
Liver_GR_WT_6pm_exo_pool/TSVrd		"mlds_605_6_CGATGT/TSVrd mlds_605_6_TGACCA/TSVrd"
Liver_GR_WT_pred_6am_exo_pool/TSVrd	"mlds_633_7_CTTGTA/TSVrd mlds_633_7_GCCAAT/TSVrd"
)

configL=(
#Liver_GR_WT_6pm_exo_pool2/TSVrd		"mlds_659_2_CGATGT/TSVrd mlds_659_2_TGACCA/TSVrd"
#Liver_GR_GRdim_6pm_exo_pool2/TSVrd	"mlds_659_2_CTTGTA/TSVrd mlds_659_2_GCCAAT/TSVrd"

Liver_GR_WT_6pm_exo_pool12/TSVrd	"Liver_GR_WT_6pm_exo_pool/TSVrd Liver_GR_WT_6pm_exo_pool2/TSVrd"
Liver_GR_GRdim_6pm_exo_pool12/TSVrd	"Liver_GR_GRdim_6pm_exo_pool/TSVrd Liver_GR_GRdim_6pm_exo_pool2/TSVrd"
)

for (( i=0;i<${#configL[@]};i=$i+2 ))
do
	desDir=${configL[$i]}
	tDirL=( ${configL[$i+1]} )
	dataDir=`dirname $desDir`
	name=$dataDir

	isDirExist ${tDirL[@]}

	echo -e "Pooling: ${desDir}" >&2
	echo ${tDirL[@]} | gawk '{ for(i=1;i<=NF;i++) printf "\t%s\n", $i }' >&2

#	continue
	mkdir -p ${desDir}
	echo -e "$name" > ${dataDir}/info.txt
	makeTagDirectory ${desDir} -d ${tDirL[@]} 2>&1 | tee ${desDir}/TSV.log
done

