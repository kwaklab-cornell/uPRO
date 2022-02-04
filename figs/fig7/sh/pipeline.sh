if false; then
# Merge pChRO data

cd bedgraph
bedtools unionbedg -i WB_pChRO_EK33_2003.mn.bedgraph  WB_pChRO_EK52_2003.mn.bedgraph  WB_pChRO_EK6_2003.mn.bedgraph  WB_pCRO_HK2_1910.mn.bedgraph \
	WB_pCRO_HK4_1910.mn.bedgraph WB_pChRO_EK45_2003.mn.bedgraph  WB_pChRO_EK53_2003.mn.bedgraph  WB_pCRO_HK1_1910.mn.bedgraph \
	WB_pCRO_HK3_1910.mn.bedgraph  WB_pCRO_SL1_1910.mn.bedgraph | \
	awk '{s=0;for(i=4;i<=NF;i++) s+=$i; print $1"\t"$2"\t"$3"\t"s;}' > WB_pChRO_all.mn.bedgraph
bedtools unionbedg -i WB_pChRO_EK33_2003.pl.bedgraph  WB_pChRO_EK52_2003.pl.bedgraph  WB_pChRO_EK6_2003.pl.bedgraph  WB_pCRO_HK2_1910.pl.bedgraph \
	WB_pCRO_HK4_1910.pl.bedgraph WB_pChRO_EK45_2003.pl.bedgraph  WB_pChRO_EK53_2003.pl.bedgraph  WB_pCRO_HK1_1910.pl.bedgraph \
	WB_pCRO_HK3_1910.pl.bedgraph  WB_pCRO_SL1_1910.pl.bedgraph | \
	awk '{s=0;for(i=4;i<=NF;i++) s+=$i; print $1"\t"$2"\t"$3"\t"s;}' > WB_pChRO_all.pl.bedgraph
cd ..

# Merged read counts
awk '$4>0{s+=$4;next}{s-=$4}END{print s}' bedgraph/WB_pChRO_all.pl.bedgraph bedgraph/WB_pChRO_all.mn.bedgraph

time ./bin/dpPRO -p bedgraph/WB_pChRO_all.pl.bedgraph -m bedgraph/WB_pChRO_all.mn.bedgraph -n net/GM0607-2.net -o  bedgraph/WB_pChRO.dpPRO.bedgraph
# 2 min 52 sec

awk '$4>0.985{print $1"\t"$2-50"\t"$3+49}' bedgraph/WB_pChRO.dpPRO.bedgraph | \
	bedtools merge -i - > bed/WB_dBTS.bed

# shared TSS
bedtools intersect -a bed/WB_dBTS.bed -b bed/combined.bed -u | wc -l

# assign promoter/enhancer/gene
./sh/dBTS_assign.sh bed/WB_dBTS.bed bed/hg38.refFlat.geneName.bed > bed/WB_assigned.bed

# cCRE analysis
./sh/mapccre.sh bed/WB_dBTS.bed bed/encodeCcre.bed > bed/WB_cCRE.bed


# deepBTS TSS quantification
./sh/dBTS_makecount.sh bed/WB_assigned.bed readcount/bglist.txt > readcount/db_rc.txt

# gene quantification
./sh/dBTS_erpkm.sh bed/hg38.refFlat.geneName.bed readcount/bglist.txt > readcount/ge_rc.txt

fi

# TFBS intersect
./sh/intersect_TFBS.sh
