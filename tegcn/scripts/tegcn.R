library(tidyr)
library(dplyr)

# Function to calculate coexpression coefficients
make_tegcn = function(db.fn = "readcount/dBTS_readcount.txt",
                      ge.fn = "readcount/gene_erpkm.txt",
                      edges.dir = "edges/",
                      cor.dir = "cor/") {
  # Load TSS and gene read counts
  db.rc = read.table(db.fn, header = T)[-1, ]
  ge.rc = read.table(ge.fn, header = T)[-1, ]
  
  # normalize reads
  normDESeq = function(x) {
    gm = apply(x, 1, function(v) exp(mean(log(v))))
    r = x/gm
    r = na.omit(r)
    sf = apply(r, 2, function(v) median(v))
    return(t(t(x)/sf))
  }
  
  db.rc[,-1] = normDESeq(db.rc[,-1])
  ge.rc[,-1] = normDESeq(ge.rc[,-1])
  
  # Load TF-enh edges
  tf.en = read.table(paste0(edges.dir, "dbts_tfbs.bed"),
                     col.names = c("chr", "start", "end", "id", "TF", "strand"))
  
  
  # Process the list of all TFs to their expression levels in gene body
  tf.list = data.frame(TF = unique(tf.en$TF))
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
  tf.list = data.frame(TF = unique(tf.en$TF))
  tf.list = tf.list %>% left_join(ge.tf %>% select(TF, replace), by = "TF") %>%
    distinct(TF, .keep_all = T) %>%
    left_join(read.table("readcount/tf_replace.txt", header = T), by = "TF") %>%
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
  tf.en = tf.en %>%
    separate(id, c("type", "gene", "geneTSS", "geneStrand"), sep = ":", remove = F) %>%
    unite("range", start, end, sep = "-", remove = F) %>%
    unite("pos", chr, range, sep = ":") %>%
    unite("id", pos, id, strand, sep = ";") %>%
    filter(TF %in% tf.list$TF) %>%
    mutate(TF = as.character(TF), id = as.character(id), gene = as.character(gene))
  
  # Calculate correlation coefficients between TSS and TF expression
  tf.en = cbind(tf.en, t(sapply(1:nrow(tf.en), function(i) {
    t = suppressWarnings(cor.test(db.mat[tf.en$id[i],], tf.mat[tf.en$TF[i],]))
    return(c(t$estimate, pval = t$p.value))})))
  tf.en = tf.en %>%
    mutate(slog.p = -sign(cor) * log10(pval))
  write.table(tf.en, paste0(cor.dir, "tf-en.txt"), row.names = F, quote = F, sep = "\t")
  
  # Calculate correlation coefficients between TSS and target gene
  # TSS to gene connection table

  db.gn = read.table(paste0(edges.dir, "enh_tss.txt"),
                     col.names = c("id", "gnid", "dist")) %>%
    separate(gnid, c("pos", "gene", "strand"), sep = ";", remove = F) %>%
    filter(gene %in% ge.tf$gene) %>%
    select(id, gene, dist) %>%
    filter(gene %in% ge.tf$gene) %>%
    filter(id %in% db.rc$id) %>%
    distinct(id, gene, .keep_all = T) %>%
    mutate(id = as.character(id), gene = as.character(gene))
  
  db.gn = cbind(db.gn, t(sapply(1:nrow(db.gn), function(i) {
    t = suppressWarnings(cor.test(db.mat[db.gn$id[i],], ge.mat[db.gn$gene[i],]))
    return(c(t$estimate, pval = t$p.value))})))
  db.gb = db.gn %>%
    mutate(slog.p = -sign(cor) * log10(pval))
  write.table(db.gn, paste0(cor.dir, "en-gn.txt"), row.names = F, quote = F, sep = "\t")
  return(
    list(db = db.mat,
         ge = ge.mat,
         tf = tf.mat,
         tf_en = tf.en,
         db_gn = db.gn)
  )
}

# Network graph library
library(network)
library(sna)
library(GGally)
library(RColorBrewer)

# enhancers with significant interactions
make_2tfnet = function(tegcn, tf1 = "E47", tf2 = "USF", chr = "chr17",
                       node.col = brewer.pal(3, "Set2"),
                       edge.col = rev(brewer.pal(11, "PiYG"))){
  tfbs.sig = tegcn$tf_en %>%
    drop_na %>%
    mutate(fdr = p.adjust(pval, method = "fdr")) %>%
    filter(fdr < 0.05)
  dbgn.sig = tegcn$db_gn %>%
    drop_na %>%
    mutate(fdr = p.adjust(pval, method = "fdr")) %>%
    filter(fdr < 0.05)
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
  edges = tibble(from = tf1, to = tf2, cor = cor(tegcn$tf[tf1,], tegcn$tf[tf2,]))
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
  names(node.col) = c("TF", "dBTS", "gene")
  nodes$color = node.col[nodes$type]
  edge.col = colorRampPalette(edge.col)(100)
  edges$color = edge.col[floor((edges$cor + 1)*50)]
  
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
  return(list(net = tfnet, palette = node.col, edgecolor = edges$color))
}

plot_tfnet = function(tfnet, seed = 1) {
  set.seed(1)
  g = ggnet2(tfnet$net,
             color = "type",
             palette = tfnet$palette,
             edge.color = tfnet$edgecolor,
             edge.alpha = 0.4,
             size = 1.5, label = "label")
  return(g)
}

# Function to find significantly associated TFs
# Non-cognate TF-dBTS co-expression as backgroun
# for each TF, count number of target dBTS

tegcn_tfcor = function(tegcn) {
  tfbs = tegcn$tf_en
  tfbs.TF.n = tfbs %>%
    group_by(TF) %>%
    summarise(count = n())
  # for each TF, randomly pick any dBTS of 1000
  tfbs.TF.rand = list()
  for(i in 1:nrow(tfbs.TF.n))
    tfbs.TF.rand[[i]] = tibble(TF = tfbs.TF.n$TF[i],
                             id = sample(unique(tfbs$id), 1000))
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
    t = suppressWarnings(cor.test(tegcn$db[tfbs.TF.rand$id[i],], tegcn$tf[tfbs.TF.rand$TF[i],]))
    return(c(t$estimate, pval = t$p.value))})))
  tfbs.rand.cor = tfbs.rand.cor %>%
    separate(id, c("pos", "name"), sep = ";", remove = F) %>%
    separate(name, c("type"), sep = ":") %>%
    mutate(type = substring(type, 1, 3))
  tfbs.cor.all = bind_rows(tfbs %>%
                             mutate(cognate = TRUE) %>%
                             mutate(type = substring(type, 1, 3)) %>%
                             select(TF, type, id, cognate, cor),
                           tfbs.rand.cor %>%
                             mutate(cognate = FALSE) %>%
                             select(TF, type, id, cognate, cor)) %>%
    drop_na
  return(tfbs.cor.all)
}

# Plot cognate vs non-cognate
plot_cognate = function(tfcor, col = brewer.pal(3, "Set2")) {
  g = ggplot(data = tfcor, aes(x = cor, col = cognate)) +
    stat_ecdf() +
    scale_color_manual(values = col) +
    theme_bw()
  return(g)
}

# Significant TFs
tfcor_sig = function(tfcor) {
  tfcor %>%
    group_by(TF, type) %>%
    summarise(mean = mean(cor[cognate]), sd = sd(cor[cognate]),
              pval = t.test(cor[cognate], cor[!cognate])$p.value) %>%
    mutate(fdr = p.adjust(pval, method = "fdr")) %>%
    filter(fdr < 0.05) %>%
    return
}

# Bimodal violin plot (guitar)
plot_guitar = function(tfcor, n = 20, type = "all", dir = "inc",
                       col.violin = brewer.pal(3, "Set2"),
                       col.hole = brewer.pal(3, "Dark2"),
                       size.hole = 3) {
  tf.sig = tfcor_sig(tfcor)
  if(dir == "inc") tf.sig = tf.sig %>% filter(mean > 0)
  else if(dir == "dec") tf.sig = tg.sig %>% filter(mean < 0)
  if(type == "prm") tf.sig = tf.sig %>% filter(type == "prm")
  else if(type == "enh") tf.sig = tf.sig %>% filter(type == "enh")
  tf.sig = tf.sig %>% arrange(-mean) %>% head(n) %>%
    select(TF, type) %>%
    unite(id, TF, type, remove = F)

  tfcor.violin = inner_join(tf.sig,
                           tfcor %>% unite(id, TF, type),
                           by = "id")
  # Process tables for ggplot
  library(mclust)
  cor_mclust_bimod = function(x) {
    m = densityMclust(x, G = 2, plot = F)
    r = c(m$parameters$mean, m$parameters$pro)
    if(r[1]>r[2]) {r[1:2] = r[2:1]; r[3:4] = r[4:3]}
    return(data.frame(bimod = c("cor1", "cor2", "pro1", "pro2"), val = r))
  }
  tfcor.bimod = tfcor.violin %>%
    group_by(TF, type) %>%
    summarise(cor_mclust_bimod(cor)) %>%
    spread(bimod, val) %>%
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
    geom_violin(aes(x = id, y = cor, fill = type), col = "grey60", size = 0.5) +
    geom_point(data = tfcor.bimod, aes(x = id, y = cor2, size = pro2, col = type)) +
    geom_point(data = tfcor.bimod, aes(x = id, y = cor2, size = pro2), pch = 21, col = "grey60") +
    scale_color_manual(values = col.hole) +
    scale_fill_manual(values = col.violin) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    guides(size = "none") +
    scale_size(range = c(0, size.hole)) +
    scale_x_discrete(labels = tf.order$TF) +
    xlab("Transcription Factor") +
    ylab("Correlation coefficient")
  return(g)
}

# Plot enhancer-gene coexpression
source("scripts/scatterPlot.R")
plot_cordist = function(tegcn) {
  dbgn.dist = tegcn$db_gn %>%
    drop_na %>%
    mutate(group = cut(dist, c(0,1000*2^(1:10),Inf))) %>%
    mutate(lower = as.numeric( sub("\\((.+),.*", "\\1", group))/1000) %>%
    mutate(group = ifelse(lower == 0, "0-2 kb",
                          ifelse(lower > 1000, "> 1 Mb",
                                 ifelse(lower > 500, "0.5-1 Mb", paste0(lower, "-", lower*2, " kb"))))) %>%
    mutate(group = factor(group, levels = c("0-2 kb", paste0(2^(0:8), "-", 2^(1:9), " kb"), "0.5-1 Mb", "> 1 Mb"))) %>%
    filter(grepl("enh", id))
  # Process tables for ggplot
  g = ggplot(data = dbgn.dist) +
    geom_hline(yintercept = 0, linetype = "dotted", size = 0.5, col = "grey60") +
    geom_violin(aes(x = group, y = cor), fill = "grey80", col = "grey60", size = 0.75, draw_quantiles = c(0.25, 0.495, 0.505, 0.75)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    xlab("Enhancer-gene distance") +
    ylab("Correlation coefficient")
  return(g)
}
