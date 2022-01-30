library(DESeq2)

# Read readcounts and use DEseq to normalize
readTable = read.table("readcount/proseq.merge.readTable.txt", header = T)
data = readTable[ , -1]
rownames(data) = make.names(readTable$id, unique = T)

metaData = data.frame(
  row.names = colnames(data),
  patient = c("P1", "P1", "P1", "P1", "P2", "P6", "P33", "P36",
              "P45", "P46", "P52", "P53", "P54", "P55"))

dds = DESeqDataSetFromMatrix(data, metaData, design = ~patient)
dds = DESeq(dds)

normDataCounts = counts(dds, norm = T)

normData = data.frame(id = rownames(normDataCounts),
                      data.frame(normDataCounts))

res = results(dds, contrast = c("patient", "P1", "P2"))

res = data.frame(id = rownames(res), res) %>%
  mutate(P2 = padj) %>%
  select(-pvalue, -padj)

patients = c("P6", "P33", "P36","P45", "P46", "P52", "P53", "P54", "P55")

for(p in patients) {
  r = results(dds, contrast = c("patient", "P1", p))
  res[,p] = r$padj
}

res = res %>%
  select(-P36)

res.sig = res %>%
  select(id, contains("P")) %>%
  filter_at(vars(contains("P")), any_vars(. < 0.05))

res.hisig = res %>%
  select(id, contains("P")) %>%
  filter_at(vars(contains("P")), any_vars(. < 0.01))

res.topsig = res %>%
  gather(indiv, pval, -id, -baseMean, -log2FoldChange, -lfcSE, -stat) %>%
  group_by(id) %>%
  summarise(min = min(pval)) %>%
  filter(min < 0.05) %>%
  ungroup %>%
  arrange(min)


########################################
# Immune related genes

immune.list = read.table("bed/immune.txt",
                         col.names = c("desc", "id")) %>%
  select(id)

immune.data = normData %>%
  inner_join(immune.list, by = "id")

immune.sig = immune.data %>%
  inner_join(res.sig %>% select(id),
             by = "id")

immune.hisig = immune.data %>%
  inner_join(res.hisig %>% select(id),
             by = "id")

all.sig = normData %>%
  inner_join(res.hisig %>% select(id),
             by = "id") %>%
  select(-P1r, -P36)

immune.norm = immune.hisig %>%
  gather(indiv, val,-id) %>%
  group_by(id) %>%
  summarise(mean= mean(val)) %>%
  ungroup() %>%
  inner_join(immune.hisig, by = "id") %>%
  mutate_at(vars(contains("P")), funs(./mean)) %>%
  select(-P1r, -P36)

min.list = c("CXCR2",
             "RSAD2",
             "SMPD3",
             "CITED2",
             "TGFBR3",
             "ILF3",
             "PSMC4",
             "ANXA1",
             "SAMSN1")

library(pheatmap)

##############################
# Figure 7B
br = 0:100/40
pdf("pdf/Fig7B.pdf", width = 3, height = 4)
pheatmap(immune.norm[,-(1:2)],
         scale = "none",
         breaks = br,
         treeheight_col = 15,
         treeheight_row = 15,
         show_rownames = F,
         border_color = NA,
         col = col)
dev.off()

####################
# Immune gene individual expression levels

IMdata = normData %>%
  filter(id %in% min.list) %>%
  mutate(id = factor(id, levels = min.list))

IMdata.long = IMdata %>%
  gather(Indiv, expression, -id) %>%
  mutate(Indiv = factor(Indiv,
                        levels = c("P1", "P1r", "P1n","P1nh",
                                   "P2", "P6", "P33", "P36",
                                   "P45", "P46", "P52", "P53",
                                   "P54", "P55")))

####################
# Figure 7
pdf("pdf/Fig7A.pdf", width = 6, height = 4)
g = ggplot(IMdata.long, aes(x = Indiv, y = expression, fill = Indiv)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), col = "grey60", size = 0.25) +
  facet_wrap( ~ id, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6)) +
  theme(strip.background = element_rect(fill="grey96")) +
  scale_fill_manual(values = col2) +
  labs(fill = "Individual",
       x = "", y = "pChRO expression")
print(g)
dev.off()


