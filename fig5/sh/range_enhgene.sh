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