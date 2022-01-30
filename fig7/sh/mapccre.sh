#!/bin/bash
# Compare a bed to cCREs and label them, Nature 583, 699–710 (2020)
# Compare HAP dpPRO to cCREs

# Parse arguments
if [ $# -lt 1 ]; then
	echo -e "Usage:\t mapccre.sh <bed> <ccre:default \"./encodeCcre.bed\">"
fi
CCRE="./encodeCcre.bed"
if [ $# -gt 1 ]; then
CCRE=$2
fi
BED=$1
TD="_tmp_mapccre"

# Make temp directory
mkdir $TD

# ccre promoters and enhancers only
awk '$4=="prom"{print $1"\t"$2"\t"$3"\t1_prom"; next}
  $4=="enhP"{print $1"\t"$2"\t"$3"\t2_enhP"; next}
  $4=="enhD"{print $1"\t"$2"\t"$3"\t3_enhD";}' $CCRE \
  > $TD/encodeEP.bed

cut -f1-3 $BED > $TD/sample.bed

# Intersect with the bed file and output
bedtools intersect -a $TD/sample.bed -b $TD/encodeEP.bed -wao \
  | cut -f1-3,7 | sort -k1,1 -k2,2n -k3,3n -k4,4 | sort -k1,1 -k2,2n -k3,3n -u \
  | awk 'substr($4,3)==""{print $1"\t"$2"\t"$3"\tnone";next}{print $1"\t"$2"\t"$3"\t"substr($4,3)}'
  
# Cleanup
rm $TD/*
rmdir $TD


