# Extract gene information from the gene clusters
awk 'NR==FNR{id=substr($4,1,length($4)-2); pos[id]=$1"\t"$2"\t"$3; strand[id]=$6;next}
	pos[$1]{print pos[$1]"\t"$1"\t"$4"\t"strand[$1]}' GencodeComprehensiveV26-hg38.bed Cluster.GencodeV26.txt > \
	Cluster.bed

sort -k1,1 -k2,2n -k3,3n Cluster.bed > tmp
mv tmp Cluster.bed

sort -k1,1 -k2,2n -k3,3n HEK293_15_segments.bed > tmp
mv tmp HEK293_15_segments.bed

# Map of chromHMM states
bedtools map -a Cluster.bed -b HEK293_15_segments.bed -c 2,3,4 -o collapse > cluster_chromHMM.txt
awk '{n=split($7,s,",");split($8,e,",");split($9,h,","); for(i=1;i<=n;++i) print $1"\t"s[i]"\t"e[i]"\t"$4"\t"$5"\t"h[i];
	if($6=="+") {print $1"\t"s[1]"\t"e[1]"\t"$4"\t"$5"\t"h[1] > "cluster_chromHMM.starts.txt";
		print $1"\t"s[n]"\t"e[n]"\t"$4"\t"$5"\t"h[n] > "cluster_chromHMM.ends.txt"; }
	else {print $1"\t"s[1]"\t"e[1]"\t"$4"\t"$5"\t"h[1] > "cluster_chromHMM.ends.txt";
		print $1"\t"s[n]"\t"e[n]"\t"$4"\t"$5"\t"h[n] > "cluster_chromHMM.starts.txt"; }}' cluster_chromHMM.txt > \
		cluster_chromHMM.all.txt
