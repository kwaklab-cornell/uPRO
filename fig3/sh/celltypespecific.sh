BED=(
"HAP1_PRO_1904.dpPRO"
"HAP1_uPRO_1910.dpPRO" 
"HEK_PRO_DMSO1_1701.dpPRO"
"HEK_PRO_DMSO2_1701.dpPRO"
"HEK_PRO_THAP1_1701.dpPRO"
"HEK_PRO_THAP2_1701.dpPRO"
"HeLa_PRO_1800.dpPRO"
"HeLa_PRO_WN_1903.dpPRO"
"HeLa_uPRO_ARF_2003.dpPRO"
"HeLa_uPRO_SSK_2003.dpPRO"
"LCL_GM18505_1500.PRO"
"LCL_GM18517_1500.PRO"
"LCL_GM18517r_1500.PRO"
"LCL_GM18520_1500.PRO"
"LCL_GM18520r_1500.PRO"
"LCL_GM18522_1500.PRO"
"LCL_GM18522r_1500.PRO"
"LCL_GM19099r_1500.PRO"
"LCL_GM19193_1500.PRO"
"LCL_GM19193r_1500.PRO"
"LCL_GM19222_1500.PRO"
"LCL_GM19222r_1500.PRO"
"LCL_GM19238_1500.PRO"
"LCL_GM19238r_1500.PRO"
"LCL_GM19239_1500.PRO"
"LCL_GM19239r_1500.PRO"
"MDM_uPRO_A_2010.dpPRO"
"MDM_uPRO_B_2010.dpPRO"
"MDM_uPRO_C_2010.dpPRO"
"PBMC_uPRO_1011.dpPRO"
"PMNL_uPRO_1911.dpPRO"
"THP1_PRO_D0_1903.dpPRO"
"THP1_PRO_D1_1903.dpPRO"
"THP1_PRO_D2_1903.dpPRO"
"THP1_PRO_D4_1903.dpPRO"
"THP1_PRO_U0_1903.dpPRO"
"THP1_PRO_U2_1903.dpPRO"
"THP1_PRO_U4_1903.dpPRO"
)
CELL=(
"HAP1"
"HAP1"
"HEK293"
"HEK293"
"HEK293"
"HEK293"
"HeLa"
"HeLa"
"HeLa"
"HeLa"
"LCL"
"LCL"
"LCL"
"LCL"
"LCL"
"LCL"
"LCL"
"LCL"
"LCL"
"LCL"
"LCL"
"LCL"
"LCL"
"LCL"
"LCL"
"LCL"
"MDM"
"MDM"
"MDM"
"PBMC"
"PMNL"
"THP1"
"THP1"
"THP1"
"THP1"
"THP1"
"THP1"
"THP1"
)

TYPE=(
"PBMC"
"PMNL"
"MDM"
"LCL"
"THP1"
"HAP1"
"BDC"
"HEK293"
"HeLa"
)

for c in "${CELL[@]}";
do
	if [ ! -f bed/${c}.cat.bed ]; then
		echo > bed/${c}.cat.bed
	fi
done

for (( i=0; i<${#BED[@]}; i++ ));
do
	cat bed/${BED[$i]}.bed >> bed/${CELL[$i]}.cat.bed
done

for c in "${CELL[@]}";
do
	if [ ! -f bed/${c}.merge.bed ]; then
		sort -k1,1 -k2,2n -k3,3n bed/${c}.cat.bed | \
			bedtools merge -i - > \
			bed/${c}.merge.bed
		rm bed/${c}.cat.bed
	fi
done

cat bed/HAP1.merge.bed \
	bed/LCL.merge.bed \
	bed/MDM.merge.bed \
	bed/PBMC.merge.bed \
	bed/PMNL.merge.bed \
	bed/THP1.merge.bed | \
	sort -k1,1 -k2,2n -k3,3n | \
	bedtools merge -i - > \
	bed/BDC.merge.bed

echo -e "chr\tstart\tend\tcCRE" > ccre/ctsTSS.txt
cut -f1-4 ccre/combined.bed >> ccre/ctsTSS.txt
for t in "${TYPE[@]}";
do
	echo ${t} > tmp
	bedtools intersect -loj -a ccre/combined.bed -b bed/${t}.merge.bed | \
		sort -k1,1 -k2,2n -k3,3n -k5,5r | \
		sort -k1,1 -k2,2n -k3,3n -u | \
		cut -f5 | \
		awk '$1=="."{print 0;next}{print 1}' >> tmp
	paste ccre/ctsTSS.txt tmp > tmp2
	mv tmp2 ccre/ctsTSS.txt
	rm tmp
done

	
