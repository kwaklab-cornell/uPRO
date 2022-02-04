t="_tmp_dir"
bed=bed/assigned.bed
tfb=bed/tfbsCons.hg38.bed

mkdir $t
sort -k1,1 -k2,2n -k3,3n $tfb > $t/tfbs.bed
bedtools intersect -wa -wb -a $bed -b $t/tfbs.bed | \
	awk '{split($10, f, /[$_]/); print $1"\t"$2"\t"$3"\t"$4"\t"f[2]"\t"$6}' | \
	sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 -k6,6 -u > \
	bed/dbts_tfbs.bed
