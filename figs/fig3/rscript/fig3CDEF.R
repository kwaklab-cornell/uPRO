ctsTSS = read.table("ccre/ctsTSS.txt", header = T)
library(FactoMineR)

ctsTSS.mca = ctsTSS %>%
  select(-chr, -start, -end, -cCRE, -BDC) %>%
  mutate_all(as.character) %>%
  MCA(ncp = 3, graph = F)
ctsTSS.var = ctsTSS.mca$var$coord %>%
  data.frame
ctsTSS.var = ctsTSS.var %>%
  mutate(cell = rownames(ctsTSS.var)) %>%
  filter(grepl("_1", cell)) %>%
  mutate(cell = substring(cell, 1, nchar(cell)-2 )) %>%
  mutate(col = ifelse(cell == "HeLa" | cell == "HEK293", "out", "blood"))

# Fig 3C, MCA analysis
pdf("pdf/Fig3C.pdf", width = 3, height = 3)
g = ggplot(data = ctsTSS.var, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(col = col), size = 2) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = col2) +
  xlab("Dim1") +
  ylab("Dim2") +
  xlim(0.25, 2.5) +
  ylim(-1, 1) +
  geom_text(label = ctsTSS.var$cell,
            nudge_y = 0.1) +
  ggtitle("deepBTS TSS - MCA")
print(g)
dev.off()

# Fig 3D: TSS overlap
ctsTSS.table = ctsTSS %>%
  group_by(BDC, HEK293, HeLa) %>%
  summarise(n = n())

bdcTSS.table = ctsTSS %>%
  mutate(type = ifelse(BDC == 1 & HEK293 == 0 & HeLa == 0, "Blood",
                       ifelse(BDC == 0, "Outgroup", "Shared"))) %>%
  group_by(type, cCRE) %>%
  summarise(n = n())

bdcTSS.frac = bdcTSS.table %>%
  ungroup %>%
  group_by(cCRE) %>%
  summarise(type = type, percent = n/sum(n)*100) %>%
  mutate(type = factor(type, levels = c("Shared", "Outgroup", "Blood")),
         cCRE = factor(cCRE, levels = c("prom", "enhP", "enhD")))

pdf("pdf/fig3D.pdf", width = 2.5, height = 3)
g = ggplot(data = bdcTSS.frac %>% filter(cCRE != "none"), aes(x = cCRE, y = percent)) +
  geom_bar(aes(fill = type), stat = "identity", position = "stack", col = "grey60") +
  scale_fill_manual(values = col[c(7,13,4)]) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
print(g)
dev.off()

# Fig S3, 3E: cell type specific TSSs within blood cells
bcTSS = ctsTSS %>%
  filter(BDC == 1) %>%
  select(-BDC, -HEK293, -HeLa) %>%
  mutate(type = PBMC+PMNL+MDM+LCL+THP1+HAP1) %>%
  mutate(type = ifelse(type >= 3, "Shared", "Specific"))

bcTSS.frac = bcTSS %>%
  group_by(type, cCRE) %>%
  summarise(n = n()) %>%
  ungroup %>%
  group_by(cCRE) %>%
  summarise(type = type, percent = n/sum(n)*100) %>%
  mutate(type = factor(type, levels = c("Shared", "Specific")),
         cCRE = factor(cCRE, levels = c("prom", "enhP", "enhD")))

bcTSS.sp = bcTSS %>%
  mutate(type = rowSums(.[5:10])) %>%
  mutate(type = ifelse(type >= 3, "Shared", "Specific")) %>%
  gather(cell, count, -chr, -start, -end, -cCRE, -type) %>%
  filter(count > 0) %>%
  select(-count) %>%
  mutate(cCRE = factor(cCRE, levels = c("none", "enhD", "enhP", "prom")))

# Fig 3E, distribution of cell type specific
pdf("pdf/fig3E.pdf", width = 4, height = 3)
g = ggplot(data = bcTSS.sp %>% filter(type == "Specific"), aes(x = cell)) +
  geom_bar(aes(fill = cCRE), position = "stack", col = "grey60") +
  scale_fill_manual(values = c("grey", col[c(13,7,4)])) +
  theme_bw() +
  ggtitle("Blood cell type specific TSSs")
print(g)
dev.off()

TSS.bc= ctsTSS%>%
  select(-BDC, -HEK293, -HeLa) %>%
  gather(cell, count, -chr, -start, -end, -cCRE) %>%
  filter(count > 0) %>%
  select(-count) %>%
  mutate(cCRE = factor(cCRE, levels = c("none", "enhD", "enhP", "prom")))

# Fig 3F, enrichment of cCREs in individual cell type specific TSSs
ctsTSS.count = ctsTSS %>%
  mutate(type = rowSums(.[5:10])) %>%
  mutate(type = ifelse(type >= 3, "Shared", "Specific")) %>%
  select(-BDC, -HEK293, -HeLa) %>%
  gather(cell, count, -chr, -start, -end, -cCRE, -type) %>%
  filter(count > 0) %>%
  select(-count) %>%
  mutate(cCRE = factor(cCRE, levels = c("none", "enhD", "enhP", "prom"))) %>%
  group_by(cell, type, cCRE) %>%
  summarise(n = n()) %>%
  ungroup %>%
  group_by(cell, cCRE) %>%
  summarise(type = type, percent = n / sum(n) * 100) %>%
  filter(type == "Specific") %>%
  ungroup %>%
  select(-type) %>%
  filter(cCRE != "none")

pdf("pdf/fig3F.pdf", width = 1.6, height = 3)
g = ggplot(data = ctsTSS.count) +
  geom_boxplot(aes(x = cCRE, y = percent, fill = cCRE), col = "grey60") +
  scale_fill_manual(values = col[c(13,7,4)]) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_discrete(limits = c("prom", "enhP", "enhD")) +
  ylab("Percent cell type specific") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
print(g)
dev.off()
