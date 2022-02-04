dbgn.dist = dbgn %>%
  mutate(group = cut(dist, c(0,1000*2^(1:10),Inf))) %>%
  mutate(lower = as.numeric( sub("\\((.+),.*", "\\1", group))/1000) %>%
  mutate(group = ifelse(lower == 0, "0-2 kb",
                        ifelse(lower > 1000, "> 1 Mb",
                               ifelse(lower > 500, "0.5-1 Mb", paste0(lower, "-", lower*2, " kb"))))) %>%
  mutate(group = factor(group, levels = c("0-2 kb", paste0(2^(0:8), "-", 2^(1:9), " kb"), "0.5-1 Mb", "> 1 Mb"))) %>%
  filter(grepl("enh", id))

# Gummyworm plot
gummyworm6 = function(tfcor.table) {
  # Process tables for ggplot
  library(mclust)
  cor_mclust_bimod = function(x) {
    m = densityMclust(x, G = 2, plot = F)
    r = c(m$parameters$mean, m$parameters$pro)
    if(r[1]>r[2]) {r[1:2] = r[2:1]; r[3:4] = r[4:3]}
    return(data.frame(bimod = c("cor1", "cor2", "pro1", "pro2"), val = r))
  }
  tfcor.bimod = tfcor.table %>%
    group_by(group) %>%
    summarise(cor_mclust_bimod(cor)) %>%
    spread(bimod, val)
  tfcor.violin = tfcor.table
  
  g = ggplot(tfcor.violin) +
    geom_hline(yintercept = 0, linetype = "dotted", size = 0.5, col = "grey60") +
    geom_violin(aes(x = group, y = cor), fill = "grey80", col = "grey60", size = 0.75, draw_quantiles = 0.5) +
    geom_point(data = tfcor.bimod, aes(x = group, y = cor2, size = pro2), col = "grey30") +
    geom_point(data = tfcor.bimod, aes(x = group, y = cor2, size = pro2), pch = 21, col = "grey60") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(size = "none") +
    xlab("Enhancer-gene distance") +
    ylab("Correlation coefficient")
  return(g)
}

g = gummyworm6(dbgn.dist)
g = g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("Correlation coefficient")


pdf("pdf/fig5F.pdf", width = 4, height = 2.5)
print(g)
dev.off()




