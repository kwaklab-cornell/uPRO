library(dplyr)
library(tidyr)

tr.s = read.table("table/sigExp.txt", header = T, stringsAsFactors = F)
pchro = read.table("table/pChRO.RPKM.txt", header = T, stringsAsFactors = F)
pchro.name = pchro %>%
	group_by(name) %>%
	summarise(Indiv1_1 = mean(HK_1),
		  Indiv1_2 = mean(HK_2),
		  Indiv1_3 = mean(HK_3),
		  Indiv2 = mean(SL_1)) %>%
        ungroup()

pchro.all = pchro.name %>%
	inner_join(tr.s, by = "name")

pchro.plot = pchro.all %>%
	mutate_at(c("Indiv1_1",
		    "Indiv1_2",
		    "Indiv2",
		    "PMNL",
		    "PBMC"),
		  funs(./pchro.all$Indiv1_3)) %>%
	select(name, Indiv1_1, Indiv1_2, Indiv2, PMNL, PBMC, sig)

library(pheatmap)
library(RColorBrewer)
col = rev(colorRampPalette(brewer.pal(11,"RdYlBu")[2:10])(101))
br = -50:50/10
pchro.sig = pchro.plot %>%
	filter(sig != "NONE") %>%
	arrange(sig)

mat = data.frame(log2(pchro.sig[,5:6]), row.names = pchro.sig$name)
annot = data.frame(Signature = factor(unlist(pchro.sig$sig)),
		   row.names = pchro.sig$name)
acol = col[c(86,15)]
names(acol) = levels(annot$Signature)
acol = list(Signature = acol)

pdf("pdf/Fig5E.pdf", width = 2.3, height = 4)
pheatmap(mat,
	cluster_rows = F,
	annotation_row = annot,
	annotation_colors = acol,
	show_rowname = F,
	cluster_cols = F,
	breaks = br,
	color = col,
	border_color = NA)
dev.off()


######################################################
pchro.avsig = pchro.plot %>%
	select(-PBMC, -PMNL) %>%
	gather(sample, val, -name, -sig) %>%
	group_by(sig, sample) %>%
	summarise(mean = 2^mean(log2(val))) %>%
	ungroup() %>%
	spread(sig, mean) %>%
	mutate(PBMC = PBMC / NONE, PMNL = PMNL / NONE) %>%
	select(-NONE) %>%
	mutate(sample = paste0("mean_", sample)) %>%
	gather(sig, mean, -sample) %>%
	spread(sample, mean) 

pchro.sesig = pchro.plot %>%
	select(-PBMC, -PMNL) %>%
	gather(sample, val, -name, -sig) %>%
	group_by(sig, sample) %>%
	summarise(se = 2^sd(log2(val)/sqrt(n()-1))) %>%
	ungroup() %>%
	spread(sig, se) %>%
	select(-NONE) %>%
	mutate(sample = paste0("se_", sample)) %>%
	gather(sig, se, -sample) %>%
	spread(sample, se) 

pchro.annot = pchro.sig %>%
	select(name, sig) %>%
	inner_join(pchro.avsig, by = "sig") %>%
	mutate_if(is.numeric, funs(5*log2(.))) %>%
	rename(Signature = sig)

col = rev(colorRampPalette(brewer.pal(11,"RdYlBu")[2:10])(101))
br = -50:50/25

mat2 = data.frame(log2(pchro.sig[,2:4]),
		 pchro.annot) %>%
	mutate(Indiv1_11 = Indiv1_1,
	       Indiv1_22 = Indiv1_2,
	       Indiv2_1 = Indiv2) %>%
	select(Indiv1_1, Indiv1_11, mean_Indiv1_1,
	       Indiv1_2, Indiv1_22, mean_Indiv1_2,
	       Indiv2, Indiv2_1, mean_Indiv2)
rownames(mat2) = pchro.annot$name

annot = data.frame(Signature = pchro.annot$Signature,
		   row.names = pchro.annot$name)

acol = col[c(86,15)]
names(acol) = levels(annot$Signature)
acol = list(Signature = acol)

pdf("pdf/Fig5F.pdf", width = 3, height = 4)
pheatmap(mat2,
	cluster_rows = F,
	annotation_row = annot,
	annotation_colors = acol,
	show_rowname = F,
	cluster_cols = F,
	breaks = br,
	color = col,
	border_color = NA,
	labels_col = c("", "Indiv1_1", "",
		       "", "Indiv1_2", "",
		       "", "Indiv2", ""),
	 gaps_col = c(3,6))
dev.off()

##################################################

refFrac = data.frame(sig = c("PBMC", "PMNL"),
		    rfr = c(0.35, 0.65))

dcp.plot = pchro.avsig %>%
	gather(sample, mean, -sig) %>%
	inner_join(refFrac, by = "sig") %>%
	mutate(fraction =  mean * rfr) %>%
	select(sig, sample, fraction) %>%
	spread(sample, fraction) %>%
	mutate_if(is.numeric, funs(./sum(.))) %>%
	gather(sample, fraction, -sig) %>%
	separate(sample, c("type", "indiv"), extra = "merge") %>%
	select(sig, indiv, fraction) %>%
	inner_join(pchro.sesig %>%
		   gather(sample, se, -sig) %>%
		   separate(sample, c("type", "indiv"), extra = "merge") %>%
		   select(sig, indiv, se),
	   by = c("sig", "indiv")) %>%
	mutate(fr_ue = se * fraction,
	       fr_le = fraction / se)

col = brewer.pal(11,"RdYlBu")[c(3,4,8)]

# Plot barplot with error bars
pdf("pdf/Fig5G.pdf", width = 3, height =3)
par(mar = c(2, 2, 1, 1), mgp = c(1.5, 0.5, 0))
ggplot(data = dcp.plot, aes(x = sig, y= fraction, fill = indiv)) +
	geom_bar(stat = "identity",
		 position = position_dodge(0.9),
		 col = "grey36",
		 size = 0.25) +
	geom_errorbar(aes(ymin = fr_le, ymax = fr_ue),
		      width = 0.4, position=position_dodge(0.9),
		      col = "black") +
	scale_fill_manual(values = col) +
	ylab("Cell type fraction") +
	xlab("Cell type") +
	coord_cartesian(ylim = c(0.25, 0.75)) +
	theme_bw() +
	theme(legend.key = element_rect())
dev.off()
