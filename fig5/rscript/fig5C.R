# Example of a TF and target enhancer
example.list = tfbs %>%
  separate(id, c("chr"), sep = ":", remove = F) %>%
  filter(TF == "E47" & type == "enhD") %>%
  arrange(-cor) %>% head(20)

k = 11

tfbs %>%
  filter(gene == example.list$gene[k] & type == "prm") 

tfbs %>%
  filter(id == example.list$id[k])

dbgn %>%
  filter(id == example.list$id[k])

scatter.example = data.frame(
  tss = db.mat[example.list$id[k],],
  gene = ge.mat[example.list$gene[k],],
  tf = tf.mat[example.list$TF[k],],
  cell = colnames(db.mat)
) %>%
  separate(cell, c("cell"), sep = "_", remove = T) %>%
  mutate(cell = factor(cell, levels = c("PBMC", "PMNL", "MDM", "LCL", "THP1")))


# Fig 5C, example scatterplot of a tf-enhancer-gene
TFname  = paste0(example.list$TF[k], "(TCF3)")
Enhname = paste0(example.list$chr[k], ":",
                 floor(example.list$start[k]/20 + example.list$end[k]/20)*10)
Genename = example.list$gene[k]
cor.example = data.frame(type = c("te", "eg", "tg"),
                         x = c(1500, 100, 1500),
                         y = rep(0, 3),
                         cor = c(cor(scatter.example$tf, scatter.example$tss),
                                 cor(scatter.example$tss, scatter.example$gene),
                                 cor(scatter.example$tf, scatter.example$gene)))
cor.example = cor.example %>%
  mutate(label = paste0("r = ", round(cor, digits = 3)))

g1 = ggplot(data = scatter.example) +
  geom_point(aes(x = tf, y = tss, col = cell)) +
  expand_limits(x = 0, y = 0) +
  theme_bw() +
  scale_color_manual(values = col4) +
  theme(legend.position = "none") +
  geom_text(data = cor.example %>% filter(type == "te"),
            aes(x = x, y = y, hjust = 0, vjust = 0,
                label = label), col = "black") +
  xlab(paste0("TF - ", TFname)) +
  ylab(paste0("Enhancer - ", Enhname))
g2 = ggplot(data = scatter.example) +
  geom_point(aes(x = tss, y = gene, col = cell)) +
  expand_limits(x = 0, y = 0) +
  theme_bw() +
  scale_color_manual(values = col4) +
  theme(legend.position = "none") +
  geom_text(data = cor.example %>% filter(type == "eg"),
            aes(x = x, y = y, hjust = 0, vjust = 0,
                label = label), col = "black") +
  xlab(paste0("Enhancer - ", Enhname)) +
  ylab(paste0("Target gene - ", Genename))
g3 = ggplot(data = scatter.example) +
  geom_point(aes(x = tf, y = gene, col = cell)) +
  expand_limits(x = 0, y = 0) +
  theme_bw() +
  geom_text(data = cor.example %>% filter(type == "tg"),
            aes(x = x, y = y, hjust = 0, vjust = 0,
                label = label), col = "black") +
  scale_color_manual(values = col4) +
  xlab(paste0("TF - ", TFname)) +
  ylab(paste0("Target gene - ", Genename))

# Generate scatterplots
library(grid)
library(gridExtra)
pdf("pdf/fig5C.pdf", width = 9, height = 2.5)
par(mar = c(2,2,2,1), mgp = c(1.5,0.5,0))
plotlist = list(g1, g2, g3)
grid.arrange(grobs=plotlist, ncol=3, widths = c(2.3, 2.4, 3.3))
dev.off()
