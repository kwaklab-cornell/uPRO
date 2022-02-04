# Fig 4C, heatmap of DE TSSs
indiv.diff = db.table %>%
  filter(CTS == T) %>%
  select(id) %>%
  inner_join(ctn.db, by = "id") %>%
  select(!ends_with("r")) %>% 
  select(!contains("HEK"))%>%
  select(!contains("HeLa"))

# Cluster count table first
clust_heatmap = function(data, k = 16, title = "") {
  data[,-1] = t(apply(data[,-1], 1, function(x) x/mean(x)))
  data = data %>% drop_na
  
  set.seed(1)
  data.clust = kmeans(data[,-1], k, iter.max = 100)
  
  data$clust = factor(paste0("c", data.clust$cluster),
                      levels = paste0("c", 1:k))
  data = data %>%
    mutate(type = ifelse(grepl("enh", id),"enh", "prm")) %>%
    arrange(type, clust) %>%
    relocate(type, clust)
  
  library(pheatmap)
  library(RColorBrewer)
  
  br = 0:100/40
  mat = as.matrix(data[,-(1:3)])
  rownames(mat) = paste0("d", 1:nrow(mat))
  acol = list(cluster = setNames(colorRampPalette(brewer.pal(11, "Paired"))(k),
                                 factor(paste0("c", 1:k), levels = paste0("c", 1:k))),
              type = setNames(col2, c("enh", "prm")))
  arow = data.frame(cluster = data$clust, type = data$type,
                    row.names = rownames(mat))
  
  pheatmap(mat, 
           breaks = br,
           color = col,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           treeheight_col = 0,
           treeheight_row = 0,
           show_rownames = F,
           show_colnames = T,
           scale = "none",
           border_color = NA,
           annotation_row = arow,
           annotation_colors = acol,
           main = title)
}

pdf("pdf/Fig4C.pdf", width = 4, height = 6)
clust_heatmap(indiv.diff, 12, title = "DE TSS (in blood cell types)")
dev.off()
