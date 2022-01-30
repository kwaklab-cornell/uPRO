# Fig5E : properties of TF-enhancer coexpression - which enhancers are better coexpressed?
# significantly co-expressed tfbs - defined in fig5A.R
tfcount.br = c(-Inf, 0, 2, 4, Inf)
tfcount.label = c("0", "1-2", "3-4", "> 4")

# Does the presence of top 20 factor binding site associated with stronger co-expression of other TFs
sigtopTF = tfbs.cor.pval %>%
  head(20) %>%
  select(TF) %>%
  unlist
sigposTF = tfbs.cor.pval %>%
  filter(mean > 0) %>%
  select(TF) %>%
  unlist
tfbs.sigcount = tfbs %>%
  group_by(id) %>%
  summarise(id = id[1], sigTFcount = sum(TF %in% sigtopTF))
tfbs.count = tfbs %>%
  filter(TF %in% sigposTF) %>%
  filter(!(TF %in% sigtopTF)) %>%
  inner_join(tfbs.sigcount, by = "id") %>%
  mutate(type = substring(type, 1,3)) %>%
  mutate(group = cut(sigTFcount, tfcount.br))

# CDF plot (for Fig S)
ggplot(tfbs.count) +
  stat_ecdf(aes(x = cor, col = group)) +
  facet_wrap(~ type)

# Gummyworm plot
gummyworm5 = function(tfcor.table, group = "group", type = "type") {
  # Process tables for ggplot
  library(mclust)
  cor_mclust_bimod = function(x) {
    m = densityMclust(x, G = 2, plot = F)
    r = c(m$parameters$mean, m$parameters$pro)
    if(r[1]>r[2]) {r[1:2] = r[2:1]; r[3:4] = r[4:3]}
    return(data.frame(bimod = c("cor1", "cor2", "pro1", "pro2"), val = r))
  }
  tfcor.bimod = tfcor.table %>%
    group_by(group, type) %>%
    summarise(cor_mclust_bimod(cor)) %>%
    spread(bimod, val) %>%
    unite(id, group, type, remove = F) %>%
    mutate(lower = as.numeric( sub("\\((.+),.*", "\\1", group) ),
          upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", group) )) %>%
    arrange(type, lower)
  tfcor.order = tfcor.bimod$id
  tfcor.bimod = tfcor.bimod %>%
    mutate(id = factor(id, levels = tfcor.order))
  tfcor.violin = tfcor.table %>%
    unite(id, group, type, remove = F) %>%
    mutate(id = factor(id, levels = tfcor.order))
  
  g = ggplot(tfcor.violin) +
    geom_hline(yintercept = 0, linetype = "dotted", size = 0.5, col = "grey60") +
    geom_violin(aes(x = id, y = cor, fill = type), col = "grey60", size = 0.75, draw_quantiles = 0.5) +
    geom_point(data = tfcor.bimod, aes(x = id, y = cor2, size = pro2, col = type)) +
    geom_point(data = tfcor.bimod, aes(x = id, y = cor2, size = pro2), pch = 21, col = "grey60") +
    # geom_point(data = tfcor.bimod, aes(x = id, y = cor1, size = pro1), pch = 21, col = "grey60") +
    scale_color_manual(values = col3) +
    scale_fill_manual(values = col2) +
    theme_bw() +
    # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    guides(size = "none") +
    xlab("Top TF binding sites") +
    ylab("Correlation coefficient") +
    scale_x_discrete(labels = rep(tfcount.label, 2))
  return(g)
}

g = gummyworm5(tfbs.count)

pdf("pdf/fig5E.pdf", width = 4, height = 2.5)
print(g)
dev.off()
