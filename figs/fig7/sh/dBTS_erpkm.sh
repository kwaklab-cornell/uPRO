#!/bin/bash
# Parse arguments
if [ $# -lt 1 ]; then
	echo -e "Usage:\t dBTS_erpkm.sh <bed> <bedgraph prefix list file>"
	exit

fi
BGLIST=$2
BED=$1
TD="_tmp_erpkm"
mkdir $TD
cut -f1-6 $BED > $TD/bed

echo id > $TD/result.txt
echo efftotal >> $TD/result.txt
awk '{print $1":"$2"-"$3";"$4";"$6}' $TD/bed >> $TD/result.txt

while read -r BG
do
	echo Processing .. $BG 1>&2
	echo $BG > $TD/col.txt
	awk '$4>0{s+=$4;c++;next}{s-=$4;c++}END{print c*c/s/1000000}' $BG.pl.bedgraph $BG.mn.bedgraph >>$TD/col.txt
	paste 	<(bedtools map -a $TD/bed -b $BG.pl.bedgraph -c 4 -o count) \
		<(bedtools map -a $TD/bed -b $BG.mn.bedgraph -c 4 -o count) \
		| awk '$6=="+"{print 0+$7;next}{print 0+$14}' >> $TD/col.txt
	paste $TD/result.txt $TD/col.txt > $TD/tmpres.txt
	mv $TD/tmpres.txt $TD/result.txt
done < $BGLIST
cat $TD/result.txt
rm -rf $TD
