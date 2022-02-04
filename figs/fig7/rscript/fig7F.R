# Gene expression table including all individuals
ge.rc = read.table("readcount/ge_rc.txt", header = T)  %>% distinct(id, .keep_all = TRUE)
samples = c("P1", "P1r", "P1n", "P1nh", "P2", "P6", "P33", "P45", "P52", "P53")
colnames(ge.rc) = c("id", samples)
pchro = ge.rc %>%
  separate(id, c("pos", "name"), sep = ";") %>%
  select(-pos) %>%
  distinct(name, .keep_all = TRUE)
pchro.name = pchro

tr.s = read.table("deconv/sigExp.txt", header = T, stringsAsFactors = F)
pchro.all = pchro.name %>%
  inner_join(tr.s, by = "name")

# Use DESeq for normalization
pchro.raw = as.matrix(pchro.all[,2:11])
rownames(pchro.raw) = pchro.all$name

coldata = data.frame(
  indiv = c(
    "P1", "P1", "P1", "P1", "P2",
    "P6", "P33", "P45", "P52", "P53"
  ),
  row.names = samples
)

library(DESeq2)
dds.raw = DESeqDataSetFromMatrix(countData = pchro.raw,
                                 colData = coldata,
                                 design = ~ indiv)
dds.raw = DESeq(dds.raw, betaPrior = T, parallel = T)

pchro.all[,2:11] = counts(dds.raw, normalized = T)
  
pchro.ref = pchro.all %>%
  select(name, P1, PBMC, PMNL) %>%
  gather(sig, expr, -name, -P1) %>%
  mutate(ref = P1) %>%
  select(-P1)

# Calculate relative transcription level between PBMC and PMNL
scale.factors = pchro.all %>%
  select(name, P1, PBMC, PMNL, sig) %>%
  group_by(sig) %>%
#  summarise(ref = mean(P1),
#            PBMC = mean(PBMC),
#            PMNL = mean(PMNL)) %>%
  summarise(ref = exp(mean(log(P1[P1>0]))),
            PBMC = exp(mean(log(PBMC[PBMC>0]))),
            PMNL = exp(mean(log(PMNL[PMNL>0])))) %>%
  gather(cell, expr, -sig, -ref)

scale.mat = scale.factors %>%
  filter(sig != "NONE") %>%
  spread(cell, expr)

rel.activity = solve(scale.mat[,3:4]) %*% as.matrix(scale.mat[,2]) %>%
  as.data.frame
rel.activity$sig = rownames(rel.activity)
rel.activity = rel.activity %>%
  inner_join(refFrac) %>%
  mutate(rel.activity = ref / rfr)

# PBMC/PMNL activity ratio
rel.activity$rel.activity[1]/rel.activity$rel.activity[2]

# Generage PBMC/PMNL fraction adjusted baseline gene expression
pchro.cal = pchro.avsig %>%
  mutate(mean_P1 = 1) %>%
  gather(indiv, fraction, -sig) %>%
  mutate(indiv = substring(indiv, 6)) %>%
  inner_join(pchro.ref, by = c("sig")) %>%
  inner_join(refFrac) %>%
  inner_join(rel.activity %>% select(sig, rel.activity)) %>%
  mutate(expr = expr * fraction * rfr * rel.activity) %>%
  select(sig, indiv, name, expr) %>%
  spread(sig, expr) %>%
  mutate(expr = PBMC + PMNL) %>%
  select(-PMNL, -PBMC) %>%
  spread(indiv, expr)

# DESeq normalization (normalize by the median of the ratio to the geometric mean)
normDESeq = function(x) {
  gm = apply(x, 1, function(v) exp(mean(log(v))))
  r = x/gm
  r = na.omit(r)
  sf = apply(r, 2, function(v) median(v))
  return(t(t(x)/sf))
}

pchro.cal = pchro.cal %>%
  relocate(P6, .after = P2) %>%
  relocate(P1r, .after = P1)

pchro.comp = pchro.all[1:11] %>%
  left_join(pchro.cal, by = "name")

pchro.norm = data.frame(gene = pchro.comp$name, normDESeq(pchro.comp[,-1]))

cor.dist = apply(pchro.norm[,-1], 1, function(x) cor(x[1:10], x[11:20]))

# Simulated cell type fraction adjusted readcounts
samplerc = read.table("readcount/ge_rc.txt", header = T, nrows = 1)[,-1] * 1000000 %>% unlist
totalcount = colSums(pchro.cal[,-1])
names(samplerc) = names(totalcount)
sample.ratio = unlist(samplerc / totalcount)
pchro.sim = round(t(t(pchro.cal[,-1]) * sample.ratio))
rownames(pchro.sim) = pchro.cal$name

# DESeq analysis
indiv.list = unique(coldata$indiv)

# Run DESeq for simulated reads
dds.sim = DESeqDataSetFromMatrix(countData = pchro.sim,
                                colData = coldata,
                                design = ~ indiv)
dds.sim = DESeq(dds.sim, betaPrior = T, parallel = T)

# Differentially expressed gene
pchro_de = function(dds, obj = "indiv") {
  diff.sim = data.frame()
  for(target in indiv.list) {
    res = results(dds,
                  contrast = list(c(paste0(obj, target)),
                                  c(paste0(obj, setdiff(indiv.list, target)))),
                  listValues = c(1, -1/(length(indiv.list) - 1))) %>%
      as.data.frame
    res$id = rownames(res)
    res = res %>%
      mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
      mutate(diff_expr = padj<0.05) %>%
      select(id, diff_expr)
    colnames(res)[2] = target
    
    if(target == indiv.list[1]) diff.sim = res
    else diff.sim = diff.sim %>% full_join(res, by = "id")
  }
  return(diff.sim)
}

pchro.de.sim = pchro_de(dds.sim)
pchro.de.raw = pchro_de(dds.raw)

# Number of DE genes in simulated data
pchro.de.sim %>%
  filter_at(vars(contains("P")), any_vars(.)) %>%
  nrow
# Number of DE genes in pChRO data
pchro.de.raw %>%
  filter_at(vars(contains("P")), any_vars(.)) %>%
  nrow

# Remove the first principal component (fractional component) of the simulated data and
# use the same eigen vector to remove the compoennt from real data
#pchro.comp[,2:11] = scale(pchro.comp[,2:11]) # scale raw data
#pchro.comp[,12:21] = scale(pchro.comp[,12:21]) # scale simulated data

# Calculate eigenvector fron the simulated data
vv = eigen(cov(pchro.comp[,12:21]))$vectors
# Rotate raw data by the eigenvector
nv = as.matrix(pchro.comp[,2:11]) %*% vv
# Remove the first PC
nv[,1] = 0
# Rotate back to data
rv = nv %*% t(vv)
# New PC1 removed,  pchro data
pchro.comp.rm = data.frame(name = pchro.comp$name, rv)
colnames(pchro.comp.rm) = colnames(pchro.cal)

pchro.pca.raw = prcomp(pchro.comp[,2:11])
pchro.pca.sim = prcomp(pchro.comp[,12:21])
pchro.pca.prm = prcomp(pchro.comp.rm[,-1])

summary(pchro.pca.raw)
summary(pchro.pca.sim)
summary(pchro.pca.prm)

sum(pchro.pca.raw$sdev^2)
sum(pchro.pca.sim$sdev^2)
sum(pchro.pca.prm$sdev^2)

pchro.pca.var = data.frame(
  sdev = c(pchro.pca.raw$sdev, pchro.pca.prm$sdev),
  PC = factor(rep(1:10, 2)),
  model = c(rep(" Raw", 10), rep("Cell fraction\ncomponent\nremoved", 10)))

#Fig 7F, Variance plot
pdf("pdf/fig7F.pdf", width = 5, height = 3)
ggplot(pchro.pca.var, aes(x = PC, y = sdev, fill = model)) +
  geom_bar(stat = "identity", position = "dodge", col = "grey60", size = 0.4) +
  theme_bw() +
  scale_fill_manual(values = col2) +
  xlab("Principal component") +
  ylab("Standard deviation")
dev.off()

