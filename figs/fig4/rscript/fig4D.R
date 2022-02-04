# PCA of all samples
ge.n.f = ge.n %>%
  filter_at(vars(-id), all_vars(. > 0.1)) %>%
  mutate_at(vars(-id), log10)

pca.ge = data.frame(prcomp(ge.n.f[,-1])$rotation[,1:2])
pca.ge$cell = rownames(pca.ge)
pca.ge = pca.ge %>%
  mutate(cell =  gsub("_.*$", "", cell)) %>%
  mutate(cell = factor(cell, levels = c("PBMC", "PMNL", "MDM", "LCL", "THP1", "HAP1", "HEK", "HeLa")))

# Fig 4E, PCA of gene body expression.
g1 = ggplot(data = pca.ge, aes(x = PC1, y = PC2)) +
  geom_point(aes(col = cell), size = 2) +
  theme_bw() +
  geom_text(label = pca.ge$cell,
            nudge_y = 0.03, check_overlap = T) +
  scale_color_manual(values = brightness(col[c(1,3,5,7,9,14,17,20)], 0.85)) +
  xlim(0.1,0.2) +
  ggtitle("Gene body expression") +
  theme(legend.position = "none")

db.en.f = db.n %>%
  filter_at(vars(-id), all_vars(. > 0.1)) %>%
  mutate_at(vars(-id), log10) %>%
  filter(grepl("enh", id))

# PCA plot in enhancers
pca.en = data.frame(prcomp(db.en.f[,-1])$rotation[,1:2])
pca.en$cell = rownames(pca.en)
pca.en = pca.en %>%
  mutate(cell =  gsub("_.*$", "", cell)) %>%
  mutate(cell = factor(cell, levels = c("PBMC", "PMNL", "MDM", "LCL", "THP1", "HAP1", "HEK", "HeLa")))

g2 = ggplot(data = pca.en, aes(x = PC1, y = PC2)) +
  geom_point(aes(col = cell), size = 2) +
  theme_bw() +
  geom_text(label = pca.ge$cell,
            nudge_y = 0.03, check_overlap = T) +
  scale_color_manual(values = brightness(col[c(1,3,5,7,9,14,17,20)], 0.85)) +
  xlim(0.1,0.2) +
  ggtitle("Enhancer TSS expression")

library(grid)
library(gridExtra)

# Fig 4D, PCA plots of enhancer and gene
pdf("pdf/Fig4D.pdf", width = 7, height = 3)
par(mar = c(2,2,2,1), mgp = c(1.5,0.5,0))
plotlist = list(g1, g2)
grid.arrange(grobs=plotlist, ncol=2, widths = c(3,4))
dev.off()

