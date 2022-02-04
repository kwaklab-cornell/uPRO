# Hierarchical TF-enhancer-gene coexpression network (TEGCN) analysis for uPRO

Nascent RNA sequencing is a powerful method to measure transcription with high resolution, sensitivity, and directional information, which gives distinctive information about transcription from other methods such as chromatin immunoprecipitation or mRNA sequencing. We present an integrated package of nascent RNA-seq methods - ultrafast Precision Run On (uPRO) combined with computational procedures to discover cell type specific enhancers, promoters, and transcription factor networks. uPRO is composed of adaptor ligation and reverse transcription reactions, which is reduced to a one-day procedure and makes nascent RNA-seq more feasible and flexible for a widespread use. We generated genome-wide profiles of nascent transcription in human blood derived cell lines and clinical samples of ~1 ml of untreated whole blood. We integrated these data into deep learning and hierarchical network analysis to detect enhancers, promoters, and co-expression networks to define cell-type specific transcription programs. We found conservation of position but variation of expression in cell type specific enhancers and transcription start sites. Transcription factors (TFs) such as TCF-3 and OCT1 were pivotally associated with TF-enhancer-gene networks across cell types. Intriguingly, we also discovered that TFs related to cell stress and inflammation - such as SRF, ATF, CHOP, and NF-kB - are associated with inter-individual variation of leukocyte transcription in whole blood. Our integration of experimental and computational nascent RNA methods will provide an efficient strategy to identify specific transcriptional programs, both in cell-type and patient/disease-associated, with minimal sample requirements.

## Running TEGCN analysis

### Required data files
- Bed file of Bidirectional Transcription Sites (BTS) positions (tegcn/bed/dBTS.bed)
- Bed file of franscription factor binding sites (TFBS) (tegcn/bed/tfbsCons.hg38.bed)
- uPRO read count table at gene bodies (tegcn/readcount/gene_erpkm.txt)
- uPRO read count table at BTS regions (tegcn/readcount/dBTS_readcount.txt)

### Contruct network edges
under UNIX shell
```
bash scripts/build_edges.sh -t <TFBS bed> -b <BTS bed> -g <gene expression table>
```
Output: network edge files under edges directory

Example
```
cd tegcn
bash ./scripts/build_edges.sh -t bed/tfbsCons.hg38.bed -b bed/dBTS.bed -g readcount/gene_erpkm.txt
```


setwd("tegcn/")
source("scripts/tegcn.R")
```

load and calculation coexpression coefficient in
```
tcn = make_tegcn()
```

Plot network on chr 17
```
tcn.chr17.net = make_2tfnet(tcn, tf1 = "E47", tf2 = "USF", chr = "chr17")
plot_tfnet(tcn.chr17.net)
```
![network plot](https://github.com/kwaklab-cornell/uPRO/blob/main/tegcn/plots/chr17_tfnet.png)


TFs-enhancer correlations that are significantly more correlated to cognate than non-cognate
```
tf.cor = tegcn_tfcor(tegcn = tcn)
tf.sig = tfcor_sig(tf.cor)
plot_guitar(tf.cor, n = 20)
```
![guitar plot](https://github.com/kwaklab-cornell/uPRO/blob/main/tegcn/plots/sigTF_guitar.png)

Distance plot
```
plot_cordist(tcn)
```
![dist plot](https://github.com/kwaklab-cornell/uPRO/blob/main/tegcn/plots/enhGene_dist.png)

### Citation
"Integrative nascent RNA methods to reveal cell-type specific transcription programs in peripheral blood and its derivative cells"

Rscripts are under figs/
