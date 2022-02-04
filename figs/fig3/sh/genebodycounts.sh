SAMPLE=(
A_M
B_M
C_M
A_RL
B_RL
C_RL
A_CD4
B_CD4
C_CD4
)
refFlat="/home/hk572/Work/shared/ref/hg38/bed/hg38.refFlat.geneName.bed"

echo id > geneBodyCounts.txt
awk '{print $1":"$2"-"$3";gb:"$4":"NR";"$6}' $refFlat >> geneBodyCounts.txt
for i in {1..9}; do
	s=${SAMPLE[$i-1]}
	echo Processing .. $s
	echo $s > _tmpColumn.txt
	paste 	<(bedtools map -a $refFlat -b bedgraph/$s.pl.bedgraph -c 4 -o count) \
		<(bedtools map -a $refFlat -b bedgraph/$s.mn.bedgraph -c 4 -o count) \
		| awk '$6=="+"{print 0+$13;next}{print 0+$26}' >> _tmpColumn.txt
	paste geneBodyCounts.txt _tmpColumn.txt > _tmpcounts.txt
	mv _tmpcounts.txt geneBodyCounts.txt
done
