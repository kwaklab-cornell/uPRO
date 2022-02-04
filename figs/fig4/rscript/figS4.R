# Coverage (eRPM) vs count (RPM)
br = 50:100/100
pdf("pdf/FigS4A.pdf", width = 6, height = 6)
pheatmap(ge.mat4, 
         treeheight_row = 20, treeheight_col = 20,
         breaks = br,
         color = col,
         main = "Gene body coverage (correlation)")
dev.off()

# Correlation of replicates
replicate.pair = list(
  c(3,4),
  c(3,5),
  c(4,5),
  c(15,16),
  c(24,25),
  c(26,27),
  c(28,29),
  c(28,30),
  c(29,30)
)

replicate.cor = data.frame(
  count = unlist(lapply(replicate.pair, function(x) gn.mat[x[1],x[2]])),
  coverage = unlist(lapply(replicate.pair, function(x) ge.mat[x[1],x[2]]))
)

pdf("pdf/FigS4B.pdf", width = 1.4, height = 3)
g = ggplot(data = replicate.cor %>% gather(type, r), aes(x = type, y = r)) +
  geom_boxplot(aes(fill = type), col = "grey60") + 
  scale_fill_manual(values = col2) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("Correlation coefficient") +
  xlab("Quantification")
print(g)
dev.off()

# PCA plot in promoters
db.pr.f = db.n %>%
  filter_at(vars(-id), all_vars(. > 0.1)) %>%
  mutate_at(vars(-id), log10) %>%
  filter(grepl("prm", id))

pca.pr = data.frame(prcomp(db.pr.f[,-1])$rotation[,1:2])
pca.pr$cell = rownames(pca.pr)
pca.pr = pca.pr %>%
  mutate(cell =  gsub("_.*$", "", cell)) %>%
  mutate(cell = factor(cell, levels = c("PBMC", "PMNL", "MDM", "LCL", "THP1", "HAP1", "HEK", "HeLa")))

# Plots in figS4C, pca in promoters
pdf("pdf/FigS4C.pdf", width = 4, height = 3)
g = ggplot(data = pca.pr, aes(x = PC1, y = PC2)) +
  geom_point(aes(col = cell), size = 2) +
  theme_bw() +
  geom_text(label = pca.pr$cell,
            nudge_y = 0.03, check_overlap = T) +
  scale_color_manual(values = brightness(col3[c(1,3,5,7,9,14,17,20)], 0.85)) +
  ggtitle("Promoter TSS expression - PCA") +
  xlim(0.125, 0.225)
print(g)
dev.off()

# Fraction of cell type variance explained
pca.var = lapply(list(ge.n.f, db.en.f, db.pr.f),
                 function(x) summary(prcomp(x[,-1]))$importance)
pca.var.name = c("Gene", "Enh", "Prm")
pca.var.table = lapply(1:3,
                       function(i) data.frame(PC = 1:8, type = pca.var.name[i], t(pca.var[[i]][,1:8])))
pca.var.table = do.call("rbind", pca.var.table) %>%
  as.data.frame

# Fig S4D, cumulative proportion in PCA
pdf("pdf/figS4D.pdf", width = 5, height = 3)
g =
ggplot(data = pca.var.table, aes(x = PC, y = Cumulative.Proportion, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", col = "grey60") +
  theme_bw() +
  scale_fill_manual(values = c(col3[9], col2)) +
  xlab("Principal components") +
  ylab("Cumulative proportion of variance")
print(g)
dev.off()
