#!/bin/bash
# Assign deepBTS bed to promoters and enhancers 
# Parse arguments
if [ $# -lt 1 ]; then
	echo -e "Usage:\t dBTS_assign.sh <deepBTS bed> <refFlat bed>"
	exit

fi
refFlat=$2
BED=$1
TD="_tmp_assign"

mkdir $TD

# Define promoter as 200 b upstream and 200 b downstream from the TSS (cCRE definition)
awk '$6=="+"{s=$2-200;if(s<1) s=1;print $1"\t"s"\t"$2+200"\t"$4"\t"$5"\t+";next}
	{s=$3-200;if(s<1) s=1;print $1"\t"s"\t"$3+200"\t"$4"\t"$5"\t-"}' \
	$refFlat | \
	sort -k1,1 -k2,2n -k3,3n -u > $TD/Promoters.bed


cut -f1-3 $BED > $TD/dBTS.bed

# Intersect promoter with dBTS bed to define active promoter
bedtools intersect -a $TD/Promoters.bed -b $TD/dBTS.bed -wa | \
	sort -k1,1 -k2,2n -k3,3n -u > $TD/activePromoters.bed

# Define active gene regions from the active promoter down to 20 kb from the poly(A) site
# Also define upstream antisense as active gene region up to 2 kb
awk 'NR==FNR{tss[$1":"$3-200":"$6]=1;next}
	$6=="+"&&tss[$1":"$2":+"]{
		s=$2-200;if(s<1) s=1;print $1"\t"s"\t"$3+20000"\t"$4"\t"$2"\t+";
		s=$2-2000;if(s<1) s=1; print $1"\t"s"\t"$2+200"\tua_"$4"\t"$2"\t-"}
	$6=="-"&&tss[$1":"$3":-"]{
		s=$2-20000;if(s<1) s=1;print $1"\t"s"\t"$3+200"\t"$4"\t"$3"\t-"
		s=$3-200;if(s<1) s=1;print $1"\t"s"\t"$3+2000"\tua_"$4"\t"$3"\t+"}' \
	$TD/activePromoters.bed $refFlat | \
	sort -k1,1 -k2,2n -k3,3n -u > $TD/activeGeneRegions.bed

# Intersect dBTS bed with active promoter to define promoter TRE on both strands
# Remove upstream antisense in bidirectional genes
bedtools intersect -a $TD/dBTS.bed -b $TD/activePromoters.bed -wa -wb | \
	sort -k1,1 -k2,2n -k3,3n | \
	awk '{ss=$1"\t"$2"\t"$3"\t1prm:"$7":"$6-200":"$9"\t0"; 
		  ua=$1"\t"$2"\t"$3"\t2prm:"$7":"$6-200":"$9"\t0"; 
		  if($9=="+") {print ss"\t+";print ua"\t-"} 
		  else {print ua"\t+";print ss"\t-"}}' | \
	sort -k1,1 -k2,2n -k3,3n -k6,6 -k4,4 | \
	sort -k1,1 -k2,2n -k3,3n -k6,6 -u | \
	awk '{$4=substr($4,2); print $1"\t"$2"\t"$3"\t1"$4"\t"$5"\t"$6}' > $TD/dBTSpromoter.bed

# Intersect dBTS with active gene regions and define intragenic or proximal enhancer TRE on the antisense strand.
bedtools intersect -a $TD/dBTS.bed -b $TD/dBTSpromoter.bed -v | \
	bedtools intersect -a - -b $TD/activeGeneRegions.bed -wa -wb | \
	sort -k1,1 -k2,2n -k3,3n -u | \
	awk '{if($9=="+") s="-"; else s="+"; 
		  if(substr($7,1,3)=="ua_") {$7=substr($7,4);et="2enhP"}
		  else if(($9=="+"&&$2<$8+2000)||($9=="-"&&$3>$8-2000)) et="3enhP";
		  else et="3enhG";
		  print $1"\t"$2"\t"$3"\t"et":"$7":"$8":"$9"\t0\t"s}' \
	> $TD/ProxIGEnhancer.bed

# Define non-overlapping dBTS as distal enhancer TRE on both strands
bedtools intersect -a $TD/dBTS.bed -b $TD/activeGeneRegions.bed -v | \
	bedtools closest -a - -b $TD/activePromoters.bed | \
	sort -k1,1 -k2,2n -k3,3n -u | \
	awk '{l=$1"\t"$2"\t"$3"\t4enhD:"$7":"$6-200":"$9"\t0"; \
		print l"\t+";print l"\t-"}' > $TD/DistEnhancer.bed

# Merge all TSSs
cat $TD/dBTSpromoter.bed $TD/ProxIGEnhancer.bed $TD/DistEnhancer.bed | \
	sort -k1,1 -k2,2n -k3,3n -k6,6 -k4,4 | \
	sort -k1,1 -k2,2n -k3,3n -k6,6 -u | \
	awk '{$4=substr($4,2); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}'
