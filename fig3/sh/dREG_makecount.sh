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

echo id > dREGcounts.txt
awk '{print $1":"$2"-"$3";"$4";"$6}' dREG.pe.bed >> dREGcounts.txt
for i in {1..9}; do
	s=${SAMPLE[$i-1]}
	echo Processing .. $s
	echo $s > _tmpColumn.txt
	paste 	<(bedtools map -a dREG.pe.bed -b bedgraph/$s.pl.bedgraph -c 4 -o sum) \
		<(bedtools map -a dREG.pe.bed -b bedgraph/$s.mn.bedgraph -c 4 -o sum) \
		| awk '$6=="+"{print 0+$7;next}{print 0-$14}' >> _tmpColumn.txt
	paste dREGcounts.txt _tmpColumn.txt > _tmpcounts.txt
	mv _tmpcounts.txt dREGcounts.txt
done
