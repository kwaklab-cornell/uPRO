library(dplyr)
library(tidyr)

dir = c("PBMC",
	"PMNL")

sample = c("PBMC", "PMNL")

file = paste0("data/", dir, "/table/expression.txt")

t.long = Reduce(bind_rows, lapply(1:2, function(i) data.frame(read.table(file[i],
					      header = T,
					      stringsAsFactors = F),
				   sample = sample[i])))
t.w = t.long %>% select(id, eRPKMgb, sample) %>%
	spread(sample, eRPKMgb)

t.c = t.long %>%
	select(id, gb, sample) %>%
	spread(sample, gb) %>%
	mutate(count = PMNL + PBMC) %>%
	select(id, count)

name = read.table(paste0("data/", dir[1], "/annotation/transcripts.bed13"),
		  header = F, stringsAsFactors = F)[,c(4, 13)]
colnames(name) = c("id", "name")
t = inner_join(name, t.w, by = "id") %>%
	inner_join(t.c, by = "id")

write.table(t, "table/leukSubtype.gb.txt", quote = F, sep = "\t", row.names = F)

tr.name = t %>%
	filter(count > 50) %>%
	group_by(name) %>%
	summarise(PMNL = mean(PMNL),
	  	PBMC = mean(PBMC),
		count = min(count)) %>%
	ungroup() %>%
	mutate(r = log2(PMNL/PBMC)) %>%
	mutate(r = r - median(r)) %>%
	arrange(-r)

library(ggplot2)
library(RColorBrewer)
col = rev(colorRampPalette(brewer.pal(11,"RdYlBu")[2:10])(101))
col = c(col[15], "#000000", col[86])

tr.plot = tr.name %>%
	mutate(sig = ifelse(r > 4, "PMNL",
			    ifelse(r < -4, "PBMC", "NONE"))) %>%
	mutate(name = factor(name, levels = name),
	       sig = factor(sig, levels = c("PMNL", "NONE", "PBMC")))

pdf("pdf/Fig5B.pdf", width = 3, height = 3)
par(mar = c(2,2,1,1), mgp = c(1.5,0.5,0))
ggplot(data = tr.plot, aes(x = name, y = r, col = sig)) +
	geom_point(stat = "identity", size = 1, stroke = 0, shape = 16) +
	theme_bw() +
	theme(axis.text.x = element_blank(),
	      axis.ticks.x = element_blank(),
	      panel.grid.major.x = element_blank(),
	      legend.position = "none") +
	scale_color_manual(values = col) +
	coord_cartesian(xlim = c(-0.02 * nrow(tr.plot), 1.02 * nrow(tr.plot)),
			ylim = c(-10, 10), expand = F)+ 
	xlab("Ordered genes") +
	ylab(expression(log[2]*PMNL/PBMC))
dev.off()

write.table(tr.plot, "table/sigExp.txt", quote = F, row.names = F, sep = "\t")
write.table(tr.plot %>% filter(sig=="PBMC") %>% select(name), "table/sigGene.pbmc.txt",
	    quote = F, row.names = F, sep = "\t")
write.table(tr.plot %>% filter(sig=="PMNL") %>% select(name), "table/sigGene.pmnl.txt",
	    quote = F, row.names = F, sep = "\t")
write.table(tr.plot %>% filter(sig=="NONE") %>% select(name), "table/sigGene.back.txt",
	    quote = F, row.names = F, sep = "\t")
###########################################################################
#

