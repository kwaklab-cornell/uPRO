#######################
# TFBS anaylsis


################################
# Gene expression quantification (Fig 4)
db.rc = read.table("readcount/db_rc.txt", header = T)
ge.rc = read.table("readcount/ge_rc.txt", header = T)  %>% distinct(id, .keep_all = TRUE)
samples = c("P1", "P1r", "P1n", "P1nh", "P2", "P6", "P33", "P45", "P52", "P53")
colnames(db.rc) = c("id", samples)
colnames(ge.rc) = c("id", samples)

# Setup DESeq to normalize
# TE count table
cts.db = db.rc[-1,-1]
rownames(cts.db) = db.rc[-1,1]
# Gene count table
cts.ge = ge.rc[-1, -1]
rownames(cts.ge) = ge.rc[-1,1]
# metadata
coldata = data.frame(
  cell = c(
    "P1", "P1", "P1", "P1", "P2",
    "P6", "P33", "P45", "P52", "P53"
  ),
  row.names = colnames(db.rc)[-1]
)
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
db.rc = data.frame(id = rownames(db.rc), db.rc)
ge.rc = counts(dds.ge, normalized = T)
ge.rc = data.frame(id = rownames(ge.rc), ge.rc)

#######################
# TFBS anaylsis
# Load TFBS on TSS (Fig 5)
tfbs = read.table("tfbs/dbts_tfbs.bed",
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
# Merge TF expression levels, exclude anything with 0
ge.tf = ge.tf %>% mutate(gene = TF) %>% select(-TF, -pos, -replace) %>% relocate(gene) %>%
  filter_if(is.numeric, all_vars(. > 0))
tf.list = tf.list %>%
  inner_join(ge.tf, by = "gene") %>%
  filter_if(is.numeric, all_vars(. > 0))
db.tf = db.rc %>%
  filter_if(is.numeric, all_vars(. > 0))

# Make readcount matrices for fast access
db.mat = as.matrix(db.tf[,-1])
rownames(db.mat) = db.tf$id
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
  filter(id %in% db.tf$id) %>%
  mutate(TF = as.character(TF), id = as.character(id), gene = as.character(gene))

# Calculate correlation coefficients between TSS and TF expression
tfbs = cbind(tfbs, t(sapply(1:nrow(tfbs), function(i) {
  t = cor.test(db.mat[tfbs$id[i],], tf.mat[tfbs$TF[i],])
  return(c(t$estimate, pval = t$p.value))})))
write.table(tfbs, "tfbs/tfbs_cor.txt", row.names = F, quote = F, sep = "\t")

# Significant co-expressions
tfbs = tfbs %>%
  drop_na %>%
  mutate(fdr = p.adjust(pval, method = "fdr")) %>%
  mutate(log.p = -sign(cor) * log10(pval))

tfbs.sig = tfbs %>%
  filter(fdr < 0.05)

# Non-cognate TF-dBTS co-expression as backgroun
# for each TF, count number of target dBTS
tfbs.TF.n = tfbs %>%
  group_by(TF) %>%
  summarise(count = n())
# for each TF, randomly pick any dBTS ~ of the count
tfbs.TF.rand = list()
for(i in 1:nrow(tfbs.TF.n))
  tfbs.TF.rand[[i]] = tibble(TF = tfbs.TF.n$TF[i],
                             id = sample(unique(tfbs$id), floor(tfbs.TF.n$count)))
tfbs.TF.rand = Reduce(bind_rows, tfbs.TF.rand)
# select non-cognate TF-dBTS pairs only
tfbs.cog = tfbs %>%
  select(TF, id) %>%
  unite(TFid, TF, id)
tfbs.TF.rand = tfbs.TF.rand %>%
  unite(TFid, TF, id, remove = F) %>%
  filter(!(TFid %in% tfbs.cog$TFid)) %>%
  select(TF, id)
# Calculate correlation coefficients
tfbs.rand.cor = cbind(tfbs.TF.rand, t(sapply(1:nrow(tfbs.TF.rand), function(i) {
  t = cor.test(db.mat[tfbs.TF.rand$id[i],], tf.mat[tfbs.TF.rand$TF[i],])
  return(c(t$estimate, pval = t$p.value))})))
tfbs.rand.cor = tfbs.rand.cor %>%
  separate(id, c("pos", "name"), sep = ";", remove = F) %>%
  separate(name, c("type"), sep = ":") %>%
  mutate(type = substring(type, 1, 3))

# p-avlue of correlation coefficient different from 0
tfbs.cor.pval = tfbs %>%
  mutate(type = substring(type, 1, 3)) %>%
  drop_na() %>%
  group_by(TF, type) %>%
  summarise(mean = mean(cor), sd = sd(cor), p = t.test(cor, tfbs.rand.cor$cor[tfbs.rand.cor$TF == TF])$p.value) %>%
  mutate(fdr = p.adjust(p, method = "fdr")) %>%
  filter(fdr < 0.05) %>%
  ungroup %>% 
  arrange(-mean)

# Significantly co-expressed TFs, top 20
tfcor.sig = tfbs %>%
  select(TF, type, cor) %>%
  mutate(type = substring(type, 1, 3)) %>%
  unite(id, TF, type) %>%
  inner_join(tfbs.cor.pval %>% head(20) %>% unite(id, TF, type) %>% select(id), by = "id") %>%
  separate(id, c("TF", "type"), remove = F)

tfbs.cor.all = bind_rows(tfcor.sig %>%
                           mutate(cognate = TRUE) %>%
                           select(TF, type, cognate, cor),
                         tfbs.rand.cor %>%
                           mutate(cognate = FALSE) %>%
                           unite(id, TF, type) %>%
                           inner_join(tfbs.cor.pval %>% head(20) %>% unite(id, TF, type) %>% select(id), by = "id") %>%
                           separate(id, c("TF", "type"), remove = F) %>%
                           select(TF, type, cognate, cor))

tfbs.rand.sig = tfbs.rand.cor %>%
  unite(id, TF, type, remove = F) %>%
  filter(id %in% unique(tfcor.sig$id))

# Estimate bimodality for all tfs
bimod_violin = function(tfcor.table) {
  # Process tables for ggplot
  library(mclust)
  cor_mclust_bimod = function(x) {
    m = densityMclust(x, G = 2, plot = F)
    r = c(m$parameters$mean, m$parameters$pro)
    if(r[1]>r[2]) {r[1:2] = r[2:1]; r[3:4] = r[4:3]}
    return(data.frame(bimod = c("cor1", "cor2", "pro1", "pro2"), val = r))
  }
  tfcor.bimod = tfcor.table %>%
    group_by(TF, type) %>%
    summarise(cor_mclust_bimod(cor)) %>%
    spread(bimod, val) %>%
    unite(id, TF, type, remove = F)
  
  tfcor.violin = tfcor.table %>%
    unite(id, TF, type, remove = F)
  tf.order = tfcor.violin %>%
    group_by(id) %>%
    summarise(TF = TF[1], type = type[1], mean = mean(cor)) %>%
    arrange(type, -mean)
  tfcor.violin = tfcor.violin %>%
    mutate(id = factor(id, levels = tf.order$id))
  tfcor.bimod = tfcor.bimod %>%
    mutate(id = factor(id, levels = tf.order$id))
  
  g = ggplot(tfcor.violin) +
    geom_hline(yintercept = 0, linetype = "dotted", size = 0.5, col = "grey60") +
    geom_violin(aes(x = id, y = cor, fill = type), col = "grey60", size = 0.75) +
    geom_point(data = tfcor.bimod, aes(x = id, y = cor2, size = pro2, col = type)) +
    geom_point(data = tfcor.bimod, aes(x = id, y = cor2, size = pro2), pch = 21, col = "grey60") +
    #    geom_point(data = tfcor.bimod, aes(x = id, y = cor1, size = pro1), pch = 21, col = "grey60") +
    scale_color_manual(values = col3) +
    scale_fill_manual(values = col2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    guides(size = "none") +
    scale_x_discrete(labels = tf.order$TF) +
    xlab("Transcription Factor") +
    ylab("Correlation coefficient")
  return(g)
}

pdf("pdf/Fig7D.pdf", width = 8, height = 2.5)
g = bimod_violin(tfcor.sig)
print(g)
dev.off()




