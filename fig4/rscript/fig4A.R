# Read readcount tables

db.rc = read.table("readcount/dBTS_readcount.txt", header = T) %>% select(-HeLa_PRO_1800) %>%
  rename_with(~gsub("_[0-9]+$", "", .x)) %>% rename_with(~gsub(".", "_", .x, fixed = T)) %>%
  rename_with(~gsub("_WN", "", .x)) %>% rename_with(~gsub("SSK", "S", .x)) %>% rename_with(~gsub("ARF", "A", .x))
gn.rc = read.table("readcount/gene_readcount.txt", header = T) %>% select(-HeLa_PRO_1800) %>% distinct(id, .keep_all = TRUE) %>%
  rename_with(~gsub("_[0-9]+$", "", .x)) %>% rename_with(~gsub(".", "_", .x, fixed = T)) %>%
  rename_with(~gsub("_WN", "", .x)) %>% rename_with(~gsub("SSK", "S", .x)) %>% rename_with(~gsub("ARF", "A", .x))
ge.rc = read.table("readcount/gene_erpkm.txt", header = T) %>% select(-HeLa_PRO_1800) %>% distinct(id, .keep_all = TRUE) %>%
  rename_with(~gsub("_[0-9]+$", "", .x)) %>% rename_with(~gsub(".", "_", .x, fixed = T)) %>%
  rename_with(~gsub("_WN", "", .x)) %>% rename_with(~gsub("SSK", "S", .x)) %>% rename_with(~gsub("ARF", "A", .x))

rpmnorm = function(rc) {
  rc[-1,-1] = t(t(rc[-1,-1])/unlist(rc[1,-1]))
  return(rc[-1,])
}
db.n = rpmnorm(db.rc)
gn.n = rpmnorm(gn.rc)
ge.n = rpmnorm(ge.rc)

db.m = db.n %>%
  select(!ends_with("r")) %>%
  filter_if(is.numeric, all_vars(. > 0))
gn.m = gn.n %>%
  select(!ends_with("r")) %>%
  filter_if(is.numeric, all_vars(. > 0))
ge.m = ge.n %>%
  select(!ends_with("r")) %>%
  filter_if(is.numeric, all_vars(. > 0))

db.mat4 = cor(scale(log10(db.m[,-1])))
gn.mat4 = cor(scale(log10(gn.m[,-1])))
ge.mat4 = cor(scale(log10(ge.m[,-1])))

library(pheatmap)
br = 50:100/100

pdf("pdf/Fig4A.pdf", width = 6, height = 6)
pheatmap(db.mat4, 
         treeheight_row = 20, treeheight_col = 20,
         breaks = br,
         color = col,
         main = "DeepBTS TSS correlation")
dev.off()

