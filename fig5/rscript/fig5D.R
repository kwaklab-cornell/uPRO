# Non-cognate TF-dBTS co-expression as backgroun
# for each TF, count number of target dBTS
tfbs.TF.n = tfbs %>%
  group_by(TF) %>%
  summarise(count = n())
# for each TF, randomly pick any dBTS 1/5 of the count
tfbs.TF.rand = list()
for(i in 1:nrow(tfbs.TF.n))
  tfbs.TF.rand[[i]] = tibble(TF = tfbs.TF.n$TF[i],
                             id = sample(unique(tfbs$id), floor(tfbs.TF.n$count/5)))
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

# Fig S5, cor in cognate is higher than non-cognate
ggplot(data = tfbs.cor.all, aes(x = cor, col = cognate)) +
  geom_density() +
  scale_color_manual(values = brightness(col2, 0.85)) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dotted", size = 0.5)

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

pdf("pdf/Fig5D.pdf", width = 8, height = 2.5)
g = bimod_violin(tfcor.sig)
print(g)
dev.off()

