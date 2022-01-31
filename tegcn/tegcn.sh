#!/bin/bash
# 

if [ $# -lt 2 ]; then
    echo -e "Usage:\t tegcn.sh [options] -g <PRO-seq filename base> -g <gene annotation bed12>"
    echo -e "Options:"
    echo -e "\t-w\tPromoter proximal range (default = 500 bp)"
    exit
fi

if [ ! -d _proseq_tmp ]; then
    mkdir _proseq_tmp
fi

RNG=500
SDIR=./
while [[ $# -ge 1 ]]; do
    key="$1"
    case $key in
        -p) PRO="$2"; shift ;;
        -g) BED="$2"; shift ;;
        -w) RNG="$2"; shift ;;
        --sdir) SDIR="$2"; shift ;;
        --default) ;;
        *) ;;
    esac
    shift
done


t="_tmp_dir"
bed=bed/assigned.bed
tfb=bed/tfbsCons.hg38.bed

mkdir $t
sort -k1,1 -k2,2n -k3,3n $tfb > $t/tfbs.bed
bedtools intersect -wa -wb -a $bed -b $t/tfbs.bed | \
	awk '{split($10, f, /[$_]/); print $1"\t"$2"\t"$3"\t"$4"\t"f[2]"\t"$6}' | \
	sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 -k6,6 -u > \
	bed/dbts_tfbs.bed
  
  t="_tmp_dir"
bed=bed/assigned.bed
ge=readcount/gene_erpkm.txt

mkdir $t

# generate bed file of active gene tsses
awk 'NR<=2{next}
  {split($1,a,";"); split(a[1],b,":");split(b[2],c,"-"); \
  if(a[3]=="+") print b[1]"\t"c[1]"\t"c[1]+1"\t"$1"\t0\t+"; \
  else print b[1]"\t"c[2]"\t"c[2]"\t"$1"\t0\t-";}' $ge | \
  sort -k1,1 -k2,2n -k3,3n -k6,6 -u > $t/genepos.bed

sort -k1,1 -k2,2n -k3,3n $bed > $t/tss.bed
bedtools window -a $t/tss.bed -b $t/genepos.bed -w 1000000 | \
awk '{pos=int(($2+$3)/2); d=pos-$8;if(d<0) d=-d; print $1":"$2"-"$3";"$4";"$6"\t"$10"\t"d}' \
  > tfbs/enh_tss.txt
