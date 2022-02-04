library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

pbmc.go = read.table("GO/GOanalysis_PBMC.txt", header = T, stringsAsFactors = F, sep = "\t")
pbmc.gs = read.table("GO/GOSlimanalysis_PBMC.txt", sep = "\t", header = T, stringsAsFactors = F)
pmnl.go = read.table("GO/GOanalysis_PMNL.txt", sep = "\t", header = T, stringsAsFactors = F)

colnames(pbmc.go)[1] = "Biological_process"
colnames(pbmc.gs)[1] = "Biological_process"
colnames(pmnl.go)[1] = "Biological_process"
pbmc = bind_rows(pbmc.go[c(1,6,7,8)], pbmc.gs[c(1,6,7,8)])
colnames(pbmc)[2:4] = c("fold_enrichment","pval", "fdr")
pmnl = pmnl.go[c(1,6,7,8)]
colnames(pmnl)[2:4] = c("fold_enrichment","pval", "fdr")

pbmc = pbmc %>% 
	group_by(Biological_process) %>%
	summarise(fold_enrichment = fold_enrichment[which.min(fdr)],
		  pval = pval[which.min(fdr)],
		  fdr = min(fdr)) %>%
	ungroup() %>% arrange(-fdr) %>%
	mutate(Biological_process = factor(Biological_process, levels = Biological_process))
pmnl = pmnl %>% arrange(-fdr) %>%
	mutate(Biological_process = factor(Biological_process, levels = Biological_process))

col = rev(colorRampPalette(brewer.pal(11,"RdYlBu")[2:10])(101))
acol = col[c(86,15)]

pdf("pdf/Fig5C1.pdf", width = 4, height = 3)
par(mar = c(2, 2, 1, 1), mgp = c(1.5, 0.4, 0))
				
ggplot(data = pbmc %>% tail(25), aes(x = Biological_process, y = -log10(pval))) +
	geom_col(fill = acol[1]) +
	ylim(0, 15) +
	coord_flip() +
	theme_bw() +
	theme(axis.text=element_text(size=6)) +
	ylab(expression(-log[10]~p)) +
	xlab("PBMC signature genes")
dev.off()

pdf("pdf/Fig5C2.pdf", width = 4, height = 3)
par(mar = c(2, 2, 1, 1), mgp = c(1.5, 0.4, 0))
col = rev(brewer.pal(11,"RdYlBu")[2:10])
ggplot(data = pmnl %>% tail(25), aes(x = Biological_process, y = -log10(pval))) +
	geom_col(fill = acol[2]) +
	ylim(0, 15) +
	coord_flip() +
	theme_bw() +
	theme(axis.text=element_text(size=6)) +
	ylab(expression(-log[10]~p)) +
	xlab("PMNL signature genes")
dev.off()
