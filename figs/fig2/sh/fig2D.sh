# Compare HEK293 dpPRO sites to dREG and chromHMM
# Compare HEK293 dpPRO to cCREs, Nature 583, 699–710 (2020)
# Compare HAP dpPRO to cCREs

cd fig2
# Locations of HEK293 dpPRO sites with respect to chromHMM and encode
awk '$4=="prom"{print $1"\t"$2"\t"$3"\t1_prom"; next}
  $4=="enhP"{print $1"\t"$2"\t"$3"\t2_enhP"; next}
  $4=="enhD"{print $1"\t"$2"\t"$3"\t3_enhD";}' bed/encodeCcre.bed \
  > loc/encodeEP.bed

bedtools intersect -a bed/HEK_PRO_DMSO1_1701.dpPRO.bed -b loc/encodeEP.bed -wao \
  | cut -f1-3,7 | sort -k1,1 -k2,2n -k3,3n -k4,4 | sort -k1,1 -k2,2n -k3,3n -u \
  | awk '{print $1"\t"$2"\t"$3"\t"substr($4,3)}' \
  > loc/HEK_dpp_enc.bed

awk '$4=="TSS"{print $1"\t"$2"\t"$3"\t1_TSS"; next}
  $4=="Enh"{print $1"\t"$2"\t"$3"\t2_Enh"; next}
  $4=="Biv"{print $1"\t"$2"\t"$3"\t3_Biv"; next}
  $4=="Tx"{print $1"\t"$2"\t"$3"\t4_Tx"; next}
  $4=="Rpt"{print $1"\t"$2"\t"$3"\t5_Rpt"; next}
  $4=="Ina"{print $1"\t"$2"\t"$3"\t6_Ina";}' bed/HEK293_red_segments.bed \
  > loc/HEK_chHMM.bed

# ChromHMM intersection
awk '{m=int(($2+$3)/2); print $1"\t"m-500"\t"m+500;}' bed/HEK_PRO_DMSO1_1701.dpPRO.bed \
  | bedtools intersect -a - -b loc/HEK_chHMM.bed -wao \
  | cut -f1-3,7 | sort -k1,1 -k2,2n -k3,3n -k4,4 | sort -k1,1 -k2,2n -k3,3n -u \
  | awk '{print $1"\t"$2"\t"$3"\t"substr($4,3)}' \
  > loc/HEK_dpp_chh.bed

awk '{m=int(($2+$3)/2); print $1"\t"m-500"\t"m+500;}' bed/HEK_dREG.bed \
  | bedtools intersect -a - -b loc/HEK_chHMM.bed -wao \
  | cut -f1-3,7 | sort -k1,1 -k2,2n -k3,3n -k4,4 | sort -k1,1 -k2,2n -k3,3n -u \
  | awk '{print $1"\t"$2"\t"$3"\t"substr($4,3)}' \
  > loc/HEK_dreg_chh.bed

# encode cCRE intersection with size restriction
awk '{m=int(($2+$3)/2); print $1"\t"m-500"\t"m+500;}' bed/HEK_PRO_DMSO1_1701.dpPRO.bed \
  | bedtools intersect -a - -b loc/encodeEP.bed -wao \
  | cut -f1-3,7 | sort -k1,1 -k2,2n -k3,3n -k4,4 | sort -k1,1 -k2,2n -k3,3n -u \
  | awk '{print $1"\t"$2"\t"$3"\t"substr($4,3)}' \
  > loc/HEK_dpp_enc2.bed

awk '{m=int(($2+$3)/2); print $1"\t"m-500"\t"m+500;}' bed/HEK_dREG.bed \
  | bedtools intersect -a - -b loc/encodeEP.bed -wao \
  | cut -f1-3,7 | sort -k1,1 -k2,2n -k3,3n -k4,4 | sort -k1,1 -k2,2n -k3,3n -u \
  | awk '{print $1"\t"$2"\t"$3"\t"substr($4,3)}' \
  > loc/HEK_dreg_enc2.bed

# Locations of HAP1 dpPRO sites with respect to encode
bedtools intersect -a bed/HAP1_PRO_1904.dpPRO.bed -b loc/encodeEP.bed -wao \
  | cut -f1-3,7 | sort -k1,1 -k2,2n -k3,3n -k4,4 | sort -k1,1 -k2,2n -k3,3n -u \
  | awk '{print $1"\t"$2"\t"$3"\t"substr($4,3)}' \
  > loc/HAP_dpp_enc.bed
  
# Locations of HAP1 dpPRO sites with respect to encode
bedtools intersect -a bed/HEK_dREG.bed -b loc/encodeEP.bed -wao \
  | cut -f1-3,9 | sort -k1,1 -k2,2n -k3,3n -k4,4 | sort -k1,1 -k2,2n -k3,3n -u \
  | awk '{print $1"\t"$2"\t"$3"\t"substr($4,3)}' \
  > loc/HEK_dreg_enc.bed
