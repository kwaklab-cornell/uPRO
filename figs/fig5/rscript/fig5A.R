# Coexpression network example

# Load TSS and gene read counts and preprocess as in fig 4.
db.rc = read.table("readcount/dBTS_readcount.txt", header = T) %>% select(-HeLa_PRO_1800) %>%
  rename_with(~gsub("_[0-9]+$", "", .x)) %>% rename_with(~gsub(".", "_", .x, fixed = T))
ge.rc = read.table("readcount/gene_erpkm.txt", header = T) %>% select(-HeLa_PRO_1800) %>% distinct(id, .keep_all = TRUE) %>%
  rename_with(~gsub("_[0-9]+$", "", .x)) %>% rename_with(~gsub(".", "_", .x, fixed = T))

# DESeq to normalize
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
db.rc = counts(dds.db, normalized = T)
db.rc = data.frame(id = rownames(db.rc), db.rc) %>% select(-contains("HAP1"), -contains("HeLa"), -contains("HEK"), -contains("PMNL"))
ge.rc = counts(dds.ge, normalized = T)
ge.rc = data.frame(id = rownames(ge.rc), ge.rc) %>% select(-contains("HAP1"), -contains("HeLa"), -contains("HEK"), -contains("PMNL"))

# Load TFBS on TSS
tfbs = read.table("bed/dbts_tfbs.bed",
                  col.names = c("chr", "start", "end", "id", "TF", "strand"))

# Process the list of all TFs to their expression levels in gene body
tf.list = data.frame(TF = unique(tfbs$TF))
# gene expression table source for TFs
ge.tf = ge.rc %>%
  separate(id, c("pos", "TF", "replace"), sep = ";")
ge.tf$mean = rowMeans(ge.tf[,-(1:3)])
# sort the table by mean expression level in the same genes, and leave only the first one (most highly expressed isoform on average)
ge.tf = ge.tf %>%
  arrange(TF, -mean) %>%
  distinct(TF, .keep_all = T) %>%
  select(-mean)

# Find TF id in ge.rc
tf.list = data.frame(TF = unique(tfbs$TF))
tf.list = tf.list %>% left_join(ge.tf %>% select(TF, replace), by = "TF") %>%
  distinct(TF, .keep_all = T) %>%
  left_join(read.table("tfbs/tf_replace.txt", header = T), by = "TF") %>%
  mutate(TF = as.character(TF), replace.y = as.character(replace.y)) %>% 
  mutate(replace.y = ifelse(!is.na(replace.x), TF, replace.y)) %>%
  mutate(gene = replace.y) %>%
  select(TF, gene)
# Merge TF expression levels
ge.tf = ge.tf %>% mutate(gene = TF) %>% select(-TF, -pos, -replace) %>% relocate(gene)
tf.list = tf.list %>%
  inner_join(ge.tf, by = "gene")

# Make readcount matrices for fast access
db.mat = as.matrix(db.rc[,-1])
rownames(db.mat) = db.rc$id
ge.mat = as.matrix(ge.tf[,-1])
rownames(ge.mat) = ge.tf$gene
tf.mat = as.matrix(tf.list[,-(1:2)])
rownames(tf.mat) = tf.list$TF

# process tfbs table to extract target gene and distance
tfbs = read.table("bed/dbts_tfbs.bed",
                  col.names = c("chr", "start", "end", "id", "TF", "strand"))
tfbs = tfbs %>%
  separate(id, c("type", "gene", "geneTSS", "geneStrand"), sep = ":", remove = F) %>%
  unite("range", start, end, sep = "-", remove = F) %>%
  unite("pos", chr, range, sep = ":") %>%
  unite("id", pos, id, strand, sep = ";") %>%
  filter(TF %in% tf.list$TF) %>%
  mutate(TF = as.character(TF), id = as.character(id), gene = as.character(gene))

# Calculate correlation coefficients between TSS and TF expression
tfbs = cbind(tfbs, t(sapply(1:nrow(tfbs), function(i) {
  t = cor.test(db.mat[tfbs$id[i],], tf.mat[tfbs$TF[i],])
  return(c(t$estimate, pval = t$p.value))})))
write.table(tfbs, "tfbs/tfbs_cor.txt", row.names = F, quote = F, sep = "\t")

# Calculate correlation coefficients between TSS and target gene
# TSS to gene connection table
# Write a function to scan all target genes within certain distance (corRange = 1000000 by default)

dbgn = read.table("tfbs/enh_tss.txt", col.names = c("id", "gnid", "dist")) %>%
  separate(gnid, c("pos", "gene", "strand"), sep = ";", remove = F) %>%
  filter(gene %in% ge.tf$gene) %>%
  select(id, gene, dist) %>%
  filter(gene %in% ge.tf$gene) %>%
  filter(id %in% db.rc$id) %>%
  distinct(id, gene, .keep_all = T) %>%
  mutate(id = as.character(id), gene = as.character(gene))

dbgn = cbind(dbgn, t(sapply(1:nrow(dbgn), function(i) {
  t = cor.test(db.mat[dbgn$id[i],], ge.mat[dbgn$gene[i],])
  return(c(t$estimate, pval = t$p.value))}))) 
write.table(dbgn, "tfbs/dbgn_cor.txt", row.names = F, quote = F, sep = "\t")

# Significant co-expressions
tfbs = tfbs %>%
  drop_na %>%
  mutate(log.p = -sign(cor) * log10(pval))
tfbs.sig = tfbs %>%
  mutate(fdr = p.adjust(pval, method = "fdr")) %>%
  filter(fdr < 0.05)
dbgn = dbgn %>%
  drop_na %>%
  mutate(log.p = -sign(cor) * log10(pval))
dbgn.sig = dbgn %>%
  mutate(fdr = p.adjust(pval, method = "fdr")) %>%
  filter(fdr < 0.05)

# Find example TF
tfcor = tfbs %>% select(TF, type, log.p)
tfcor.test = tfcor %>%
  group_by(TF, type) %>%
  summarise(mean = mean(log.p), p = t.test(log.p)$p.value) %>%
  arrange(type, p)

# Test cooperativity of 2 TFs with same TSS shared TFs
test_coop = function(tf1, tf2, type ="all") {
  tf1tss = tfbs %>%
    filter(TF == tf1) %>%
    select(id) %>%
    unlist %>%
    as.character
  tf2tss = tfbs %>%
    filter(TF == tf2) %>%
    select(id) %>%
    unlist %>%
    as.character
  common.tss = intersect(tf1tss, tf2tss)
  if(type != "all") common.tss = as.character(tfbs$id[tfbs$id %in% common.tss & tfbs$type == type])
  tf1.exp = scale(tf.mat[tf1,])
  tf2.exp = scale(tf.mat[tf2,])
  tss.exp = scale(db.mat[common.tss,])
  return(data.frame(tf1 = tf1, tf2 = tf2, t(apply(tss.exp, 1, function(x) {
    data = data.frame(tf1 = tf1.exp, tf2 = tf2.exp,
               tss = x)
    m = lm(tss ~ tf1 + tf2 + tf1*tf2, data)
    return(c(m$coefficients[-1], car::vif(m)))
  }))))
}

tftestlist = tfcor.test %>%
  arrange(p) %>%
  select(TF) %>%
  unlist %>%
  unique

coop.res = list()
for(i in 1:50) if(tftestlist[i] != "E47") coop.res[[i]] = test_coop("E47", tftestlist[i])
coop.res = Reduce(bind_rows, coop.res)
coop.res.sig = coop.res %>%
  filter(tf1.1 < 10 & tf2.1 < 10 & tf1.tf2.1 < 10)

# Network graph library
library(network)
library(sna)
library(GGally)

# enhancers with significant interactions, shared between E47 and USF
make_2tfnet = function(tf1 = "E47", tf2 = "USF", chr = "chr17") {
enh.sig.id1 = tfbs.sig %>%
  filter(TF == tf1) %>%
  select(id) %>%
  filter(grepl(chr, id)) %>%
  unlist %>%
  as.character
enh.sig.id2 = tfbs.sig %>%
  filter(TF == tf2) %>%
  select(id) %>%
  filter(grepl("chr17", id)) %>%
  unlist %>%
  as.character
enh.common = intersect(enh.sig.id1, enh.sig.id2)
enh.sig.id1 = setdiff(enh.sig.id1, enh.common)
enh.sig.id2 = setdiff(enh.sig.id2, enh.common)
enh.sig.id1 = tibble(id = enh.sig.id1)
enh.sig.id2 = tibble(id = enh.sig.id2)
enh.common = tibble(id = enh.common)
enh.sig.id1 = enh.sig.id1 %>%
  inner_join(tfbs.sig %>% filter(TF == tf1) %>% select(id, cor))
enh.sig.id2 = enh.sig.id2 %>%
  inner_join(tfbs.sig %>% filter(TF == tf2) %>% select(id, cor))
enh.common = enh.common %>%
  inner_join(tfbs.sig %>% filter(TF == tf1) %>% select(id, cor)) %>%
  inner_join(tfbs.sig %>% filter(TF == tf2) %>% select(id, cor), by = "id")
dbgn.net = dbgn.sig %>%
  filter(id %in% c(enh.sig.id1$id,
                   enh.sig.id2$id,
                   enh.common$id)) %>%
  select(id, gene, cor)


# constructing network for visualization
# TFs
nodes = tibble(id = c(tf1, tf2), type = "TF")
edges = tibble(from = tf1, to = tf2, cor = cor(tf.mat[tf1,], tf.mat[tf2,]))
# Enhancer and promoters to TF1 in chr17
nodes = nodes %>% 
  bind_rows(tibble(id = enh.sig.id1$id, type = "dBTS"))
edges = edges %>%
  bind_rows(tibble(from = tf1, to = enh.sig.id1$id,
                   cor = enh.sig.id1$cor))
# Enhancer and promoters to TF2 in chr17
nodes = nodes %>% 
  bind_rows(tibble(id = enh.sig.id2$id, type = "dBTS"))
edges = edges %>%
  bind_rows(tibble(from = tf2, to = enh.sig.id2$id,
                   cor = enh.sig.id2$cor))
# Enhancer and promoters to both TFs in chr17
nodes = nodes %>% 
  bind_rows(tibble(id = enh.common$id, type = "dBTS"))
edges = edges %>%
  bind_rows(tibble(from = tf1, to = enh.common$id,
                   cor = enh.common$cor.x)) %>%
  bind_rows(tibble(from = tf2, to = enh.common$id,
                   cor = enh.common$cor.y))
# Add gene body expression
nodes = nodes %>%
  bind_rows(tibble(id = as.character(unique(dbgn.net$gene)), type = "gene"))
edges = edges %>%
  bind_rows(tibble(from = dbgn.net$id, to = dbgn.net$gene,
                   cor = dbgn.net$cor))

# Generate colors for the nodes and edges
node.color = c("#000000", col2[2], col2[1])
names(node.color) = c("TF", "dBTS", "gene")
nodes$color = node.color[nodes$type]
edges$color = col[floor((edges$cor + 1)*50)]

# fix node id
node.id = 1:nrow(nodes)
names(node.id) = nodes$id
nodes  = nodes %>%
  mutate(label = id) %>%
  mutate(id = node.id)
edges = edges %>%
  mutate(to = node.id[to], from = node.id[from])
nodes = nodes %>%
  mutate(label = ifelse(type == "TF", paste0("\t        ", label), NA))
# Make network
tfnet = network(edges,
                vertex.attr = nodes,
                matrix.type = "edgelist")
return(list(net = tfnet, palette = node.color, edgecolor = edges$color))
}

tfnet = make_2tfnet(tf1 = "E47", tf2 = "OCT")

# plot
pdf("pdf/Fig5A.pdf", width = 5, height = 4)
set.seed(1)
g = ggnet2(tfnet$net,
       color = "type",
       palette = tfnet$palette,
       edge.color = tfnet$edgecolor,
       edge.alpha = 0.4,
      size = 1.5, label = "label")
print(g)
dev.off()

