# DESeq to identify differentially expressed TEs and genes
# TE count table
cts.db = db.rc[-1,-1]
rownames(cts.db) = db.rc[-1,1]
# Gene count table
cts.ge = ge.rc[-1, -1]
rownames(cts.ge) = ge.rc[-1,1]

# Cell type data
coldata = data.frame(
  cell = c(
    rep("PBMC", 1),
    rep("PMNL", 1),
    rep("MDM", 3),
    rep("LCL", 17),
    rep("HAP1", 2),
    rep("THP1", 7),
    rep("HEK293", 4),
    rep("HeLa", 3)
  ),
  row.names = colnames(db.rc)[-1]
)
cell.list = unique(as.character(coldata$cell))

# Run DESeq
library(DESeq2)
dds.db = DESeqDataSetFromMatrix(countData = cts.db,
                                colData = coldata,
                                design = ~ cell)
dds.db = DESeq(dds.db, betaPrior = T, parallel = T)
dds.ge = DESeqDataSetFromMatrix(countData = cts.ge,
                                colData = coldata,
                                design = ~ cell)
dds.ge = DESeq(dds.ge, betaPrior = T, parallel = T)

# Normalized counts
ctn.db = counts(dds.db, normalized = T)
ctn.db = data.frame(id = rownames(ctn.db), ctn.db)
ctn.ge = counts(dds.ge, normalized = T)
ctn.ge = data.frame(id = rownames(ctn.ge), ctn.ge)

# Blood derived cells vs the outgroup
cell.bd = cell.list[1:6]
cell.og = cell.list[7:8]

compare.outgroup = function(dds, group, outgroup) {
  ct.table = data.frame()
  for(target in group) {
    res = results(dds,
                  contrast = list(c(paste0("cell", target)),
                                  c(paste0("cell", outgroup))),
                  listValues = c(1, -1/(length(outgroup)))) %>% as.data.frame
    res$id = rownames(res)
    res = res %>%
      mutate(cell_specific = ifelse(abs(log2FoldChange)>1 & padj<0.05, T, F)) %>%
      select(id, cell_specific)
    colnames(res)[2] = target
    if(target == group[1]) ct.table = res
    else ct.table = ct.table %>% full_join(res, by = "id")
  }
  ct.table$UNIT = apply(ct.table[,-1], 1, any)
  res = results(dds,
                contrast = list(c(paste0("cell", group)),
                                c(paste0("cell", outgroup))),
                listValues = c(1/length(group), -1/(length(outgroup)))) %>% as.data.frame
  res$id = rownames(res)
  res = res %>%
    mutate(MEAN = ifelse(abs(log2FoldChange)>1 & padj<0.05, T, F)) %>%
    select(id, MEAN)
  ct.table = ct.table%>% full_join(res, by = "id")
  ct.table[is.na(ct.table)] = F
  return(ct.table)
}

# Blood derived cell type specific TSSs and gene expression
db.ct.bd = compare.outgroup(dds.db, cell.bd, cell.og)
ge.ct.bd = compare.outgroup(dds.ge, cell.bd, cell.og)

# Individual cell type specific TSSs and gene expression
compare.ingroup = function(dds, group) {
  ct.table = data.frame()
  for(target in group) {
    res = results(dds,
                  contrast = list(c(paste0("cell", target)),
                                  c(paste0("cell", setdiff(group, target)))),
                  listValues = c(1, -1/(length(group)-1))) %>% as.data.frame
    res$id = rownames(res)
    res = res %>%
      mutate(cell_specific = ifelse(abs(log2FoldChange)>1 & padj<0.05, T, F)) %>%
      select(id, cell_specific)
    colnames(res)[2] = target
    if(target == group[1]) ct.table = res
    else ct.table = ct.table %>% full_join(res, by = "id")
  }
  ct.table[is.na(ct.table)] = F
  ct.table$CTS = apply(ct.table[,-1], 1, any)
  ct.table$COM = !ct.table$CTS
  return(ct.table)
}

# Cell type specificity table
ge.table = compare.ingroup(dds.ge, cell.bd)
db.table = compare.ingroup(dds.db, cell.bd) 

# Some statistics
# Fraction of differentially expressed genes
diff_frac = function(table, CELL) {
  df = data.frame(type = c("all" ,"enhD", "enhP", "enhG", "prm"),
                  frac = c(table %>% filter(!!sym(CELL)==T) %>% nrow /
                             table %>% nrow * 100,
                           table %>% filter(!!sym(CELL)==T & grepl("enhD", id)) %>% nrow /
                             table %>% filter(grepl("enhD", id)) %>% nrow * 100,
                           table %>% filter(!!sym(CELL)==T & grepl("enhP", id)) %>% nrow /
                             table %>% filter(grepl("enhP", id)) %>% nrow * 100,
                           table %>% filter(!!sym(CELL)==T & grepl("enhG", id)) %>% nrow /
                             table %>% filter(grepl("enhG", id)) %>% nrow * 100,
                           table %>% filter(!!sym(CELL)==T & grepl("prm", id)) %>% nrow /
                             table %>% filter(grepl("prm", id)) %>% nrow * 100)) %>% drop_na
  return(df)
}

# Conver db.table to a long list
db.long = db.table %>%
  gather(cell, specific, -id) %>%
  separate(id, c("pos", "id", "strand"), sep = ";") %>%
  separate(id, c("type", "gene", "center", "strand"), sep = ":") %>%
  mutate(type = substring(type, 1, 3))%>%
  filter(cell != "COM" & cell != "CTS")

db.t = db.long %>%
  filter(specific == T) %>%
  distinct(pos, cell, gene, .keep_all = T)

pdf("pdf/Fig4B.pdf", width = 3, height = 3)
g = ggplot(data = db.t, aes(x = cell)) +
  geom_bar(aes(fill = type), position = "stack", col = "grey60") +
  scale_fill_manual(values = col2) +
  theme_bw() +
  labs(fill = "TSS type") +
  ggtitle("DE TSSs") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
print(g)
dev.off()
