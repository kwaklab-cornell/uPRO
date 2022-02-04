# Cell type specific gene expression analysis
library(dplyr)
library(tidyr)

samples = c("HAP1", "HEK.rep1", "HEK.rep2",
            "HeLa", "nuHeLa")
geneList = read.table("allCells/expr/gencode.v26.geneName.txt")
geneNames = geneList[,c(4,13)]
geneLen = data.frame(id = geneList[,4], len = (geneList[,3] - geneList[,2])/1000)

tables = lapply( paste0("allCells/expr/", samples, ".expr.txt"), read.table, header = 1)
names(tables) = samples

eRPKMgb = data.frame(id = geneNames[,1], geneName = geneNames[,2], sapply(tables, function(x) return(x$eRPKMgb)))
RPKMpp = data.frame(id = geneNames[,1], geneName = geneNames[,2], sapply(tables, function(x) return(x$RPKMpp)))

LCL = read.table("allCells/expr/LCLgb.txt", header = T, stringsAsFactors = F)
LCLtrc = c(22.29, 23.76, 16.51, 24.45)
LCL[,-1] = t(t(LCL[,-1])/LCLtrc)
LCLgb = inner_join(geneLen, LCL, by = "id") %>%
       mutate_at(vars(contains("GM")), funs(./len)) %>%
       select(-len)

eRPKMgb2 = eRPKMgb %>%
	inner_join(LCLgb, by = "id") %>%
	select(id, HEK.rep1, HEK.rep2, HeLa, nuHeLa, GM18520.r2, GM19222.r2)

eRPKMgb %>%
  filter(HAP1 > 1 & HEK.rep1 < 0.1 & HEK.rep2 < 0.1 &
           HeLa < 0.1) %>%
  select(geneName) %>%
  unique() %>%
  left_join(eRPKMgb) ->
  HAP1_specific_genes


# PRO-seq vw nuPRO correlation
samples = c("PRO_1904", "nuPRO_1910")
files = paste0("HAP1/", samples, "/table/expression.txt")

data.raw = Reduce(bind_rows,
		  lapply(1:2,
			 function(i) data.frame(read.table(files[i],
							   header = T,
							   stringsAsFactors = F),
						sample = samples[i])))
data.gb = data.raw %>%
	select(id, eRPKMgb, sample) %>%
	separate(sample, c("sample")) %>%
	spread(sample, eRPKMgb) %>%
	filter_if(is.numeric, all_vars(. > 0))
data.pp = data.raw %>%
	select(id, RPKMpp, sample) %>%
	separate(sample, c("sample")) %>%
	spread(sample, RPKMpp) %>%
	filter_if(is.numeric, all_vars(. > 0))

source("rscript/scatterPlot.R")
library(ggplot2)
library(grid)
library(gridExtra)

pdf("pdf/Fig1B.pdf", width = 6, height = 3)
par(mar = c(2,2,2,1), mgp = c(1.5,0.5,0))
plotlist = list(plot.cor(data.frame(x = data.pp[,2], y = data.pp[,3]),
						"uPRO(HAP1)", "PRO-seq(HAP1)", cor=T, diag=F, red=100) +
		labs(title = "Promoter proximal"),
		plot.cor(data.frame(x = data.gb[,2], y = data.gb[,3]),
						"uPRO(HAP1)", "PRO-seq(HAP1)", cor=T, diag=F, red=100) +
		labs(title = "Gene body"))
grid.arrange(grobs=plotlist, ncol=2)
dev.off()

pdf("pdf/Fig1C.pdf", width = 6, height =3)
par(mar = c(2,2,2,1), mgp = c(1.5,0.5,0))
plotlist = list(plot.cor(data.frame(x = RPKMpp$nuHeLa, y = RPKMpp$HeLa),
						"uPRO(HeLa)", "PRO-seq(HeLa)", cor=T, diag=F, red=100) +
		labs(title = "Promoter proximal"),
		plot.cor(data.frame(x = eRPKMgb$nuHeLa, y = eRPKMgb$HeLa),
						"uPRO(HeLa)", "PRO-seq(HeLa)", cor=T, diag=F, red=100) +
		labs(title = "Gene body"))
grid.arrange(grobs=plotlist, ncol=2)
dev.off()

pdf("pdf/Fig1D.pdf", width = 6, height =3)
par(mar = c(2,2,2,1), mgp = c(1.5,0.5,0))
plotlist = list(plot.cor(data.frame(x = RPKMpp$HEK.rep1, y = RPKMpp$HEK.rep2),
						"HEK293(rep1)", "HEK293(rep2)", cor=T, diag=F, red=100) +
		labs(title = "Promoter proximal"),
		plot.cor(data.frame(x = eRPKMgb$HEK.rep1, y = eRPKMgb$HEK.rep2),
						"HEK293(rep1)", "HEK293(rep2)", cor=T, diag=F, red=100) +
		labs(title = "Gene body"))
grid.arrange(grobs=plotlist, ncol=2)
dev.off()

gb.all = eRPKMgb2 %>%
	inner_join(data.gb, by = "id") %>%
	mutate(HAP1.nuPRO = nuPRO, HAP1.PRO = PRO) %>%
	select(-nuPRO, -PRO) %>%
	filter_if(is.numeric, all_vars(. > 0))
mat = cor(scale(log10(gb.all[,-1])))

library(pheatmap)
library(RColorBrewer)

br = 50:100/100
pdf("pdf/Fig1E.pdf", width = 3.5, height = 3)
pheatmap(mat, labels_row = c("HEK293(rep1)",
	"HEK293(rep2)",
	"HeLa(PRO-seq)",
	"HeLa(uPRO)",
	"GM18520",
	"GM19222",
	"HAP1(uPRO)",
	"HAP1(PRO-seq)"),
	labels_col = c("HEK.1", "HEK.2", "HeLa.PRO","HeLa.uP","GM18520","GM19222","HAP1.uP","HAP1.PRO"),
	treeheight_row = 20, treeheight_col = 20,
	breaks = br,
	color = col,
	main = "Gene body correlation")
dev.off()
