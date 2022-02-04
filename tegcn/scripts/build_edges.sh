#!/bin/bash
# Build TEGCN TF-enhancer and enhancer-gene edges.

# Print usage
if [ $# -lt 2 ]; then
  echo -e "Usage:\t build_edges.sh [options] -t <TFBS bed> -b <BTS bed> -g <gene expression table>"
	echo -e "Options:"
	echo -e "\t--hr\tHigh resolution version"
	echo -e "\t--out\tOutput directory (default = \"edges\")"
	exit
fi

# Temp directory
if [ ! -d _tmp_dir ]; then
  mkdir _tmp_dir
fi
TD="_tmp_dir"

# Parse arguments
HR=false
OUT="edges"
while [[ $# -ge 1 ]]; do
	key="$1"
	case $key in
		-t) TF="$2"; shift ;;
		-b) BTS="$2"; shift ;;
		-g) GE="$2"; shift ;;
		--hr) HR=true; ;;
		--out) OUT="$2"; shift ;;
		--default) ;;
		*) ;;
	esac
	shift
done

# Build TF-BTS edges
sort -k1,1 -k2,2n -k3,3n $TF > $TD/tfbs.bed
bedtools intersect -wa -wb -a $BTS -b $TD/tfbs.bed | \
	awk '{split($10, f, /[$_]/); print $1"\t"$2"\t"$3"\t"$4"\t"f[2]"\t"$6}' | \
	sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 -k6,6 -u > \
	$OUT/dbts_tfbs.bed

# generate position bed file of active gene tss's
awk 'NR<=2{next}
  {split($1,a,";"); split(a[1],b,":");split(b[2],c,"-"); \
  if(a[3]=="+") print b[1]"\t"c[1]"\t"c[1]+1"\t"$1"\t0\t+"; \
	else print b[1]"\t"c[2]"\t"c[2]"\t"$1"\t0\t-";}' $GE | \
	sort -k1,1 -k2,2n -k3,3n -k6,6 -u > $TD/genepos.bed

# Build all edges between BTS and genes within 1 Mb
sort -k1,1 -k2,2n -k3,3n $BTS > $TD/tss.bed
bedtools window -a $TD/tss.bed -b $TD/genepos.bed -w 1000000 | \
awk '{pos=int(($2+$3)/2); d=pos-$8;if(d<0) d=-d; print $1":"$2"-"$3";"$4";"$6"\t"$10"\t"d}' \
  > $OUT/enh_tss.txt

rm $TD/*
rmdir $TD
