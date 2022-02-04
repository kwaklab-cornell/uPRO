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

# Intersect promoter with dBTS bed to define active promoter
bedtools intersect -a $TD/Promoters.bed -b $BED -wa | \
	sort -k1,1 -k2,2n -k3,3n -u > $TD/activePromoters.bed

# Define active gene regions from the active promoter down to 20 kb from the poly(A) site
# Also define upstream antisense as active gene region up to 2 kb
awk 'NR==FNR{tss[$1":"$3-200":"$6]=1;next}
	$6=="+"&&tss[$1":"$2":+"]
		{s=$2-200;if(s<1) s=1;print $1"\t"s"\t"$3+20000"\t"$4"\t"$5"\t+";
		 s=$2-2000;if(s<1) s=1; print $1"\t"s"\t"$2+200"\tua_"$4"\t"$5"\t-";next}
	$6=="-"&&tss[$1":"$3":-"]
		{s=$2-20000;if(s<1) s=1;print $1"\t"s"\t"$3+200"\t"$4"\t"$5"\t-"
		 s=$3-200;if(s<1) s=1;print $1"\t"s"\t"$3+2000"\tua_"$4"\t"$5"\t+"}' \
	$TD/activePromoters.bed $refFlat | \
	sort -k1,1 -k2,2n -k3,3n -u > $TD/activeGeneRegions.bed

# Intersect dBTS bed with active promoter to define promoter TRE on both strands
bedtools intersect -a $BED -b $TD/activePromoters.bed -wa -wb | \
	sort -k1,1 -k2,2n -k3,3n -u | \
	awk '{l=$1"\t"$2"\t"$3"\tprm:"$9":"$11"\t"$5; \
		print l"\t+";print l"\t-"}' > $TD/dBTSpromoter.bed

# Intersect dBTS with active gene regions and define intragenic or proximal enhancer TRE on the antisense strand.
bedtools intersect -a $BED -b $TD/dBTSpromoter.bed -v | \
	bedtools intersect -a - -b $TD/activeGeneRegions.bed -wa -wb | \
	sort -k1,1 -k2,2n -k3,3n -u | \
	awk '{if($11=="+") s="-"; else s="+"; 
		  if(substr($9,1,3)=="ua_") {$9=substr($9,4);et="enhP"} else et="enhG";
		  print $1"\t"$2"\t"$3"\t"et":"$9":"$11"\t"$5"\t"s}' \
	> $TD/ProxIGEnhancer.bed

# Define non-overlapping dBTS as distal enhancer TRE on both strands
bedtools intersect -a $BED -b $TD/activeGeneRegions.bed -v | \
	bedtools closest -a - -b $TD/activePromoters.bed | \
	sort -k1,1 -k2,2n -k3,3n -u | \
	awk '{l=$1"\t"$2"\t"$3"\tenhD:"$9":"$11"\t"$5; \
		print l"\t+";print l"\t-"}' > $TD/DistEnhancer.bed

# Merge all dREGs
cat $TD/dBTSpromoter.bed $TD/ProxIGEnhancer.bed $TD/DistEnhancer.bed | \
	sort -k1,1 -k2,2n -k3,3n -k6,6
