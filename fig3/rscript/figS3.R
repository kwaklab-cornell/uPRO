# Fig S3

# Fig S3A, number of each cCREelements
ccre.table = ccre.table %>%
  sub("SSK", "S", sample) %>%
  sub("ARF", "A", sample)
pdf("pdf/FigS3A.pdf", width = 4, height = 6)
g = ggplot(data = ccre.table, aes(x = sample)) +
  geom_bar(aes(fill = cCRE), position = "stack", col = "grey60", size = 0.1) +
  scale_fill_manual(values = c("grey", col[c(13,7,4)])) +
  coord_flip() + 
  theme_bw()
print(g)
dev.off()

pdf("pdf/FigS3B.pdf", width = 2.5, height = 3)
g = ggplot(data = bcTSS.frac %>% filter(cCRE != "none"), aes(x = cCRE, y = percent)) +
  geom_bar(aes(fill = type), stat = "identity", position = "stack", col = "grey60") +
  scale_fill_manual(values = col[c(7,4)]) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
print(g)
dev.off()

pdf("pdf/FigS3C.pdf", width = 3.5, height = 3)
g = ggplot(data = TSS.bc, aes(x = cell)) +
  geom_bar(aes(fill = cCRE), position = "stack", col = "grey60") +
  scale_fill_manual(values = c("grey", col[c(13,7,4)])) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle("All TSSs")
print(g)
dev.off()
