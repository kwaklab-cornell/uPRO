awk '{printf $1"\t"$2"\t"$3"\t"}
	$4=="E1"||$4=="E2"{print "TSS";next}
	$4=="E3"||$4=="E4"||$4=="E5"{print "Tx";next}
	$4=="E6"||$4=="E7"{print "Enh";next}
	$4=="E8"{print "Rpt"; next}
	$4=="E9"||$4=="E13"||$4=="E14"||$4=="E15"{print "Ina"; next}
	$4=="E10"||$4=="E11"||$4=="E12"{print "Biv"; next}' HEK293_15_segments.bed | \
	awk 'NR==1{prevChr=$1; prevStart=$2; prevEnd=$3; prevState=$4;}
		$4!=prevState||$1!=prevChr{
			print prevChr"\t"prevStart"\t"prevEnd"\t"prevState;
			prevChr=$1; prevStart=$2; prevEnd=$3; prevState=$4;
			next;}
		{	prevEnd=$3;	}'


