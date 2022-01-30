# Read HEK and HAP dBTS regions maps
dbts.map = list()

dbts.map$HAP = read.table("loc/HAP_dpp_enc.bed", fill = T,
                          col.names = c("chr", "start", "end", "cCRE")) %>%
  mutate(cCRE = if_else(cCRE=="", "none", as.character(cCRE))) 
dbts.map$HEK_enc = read.table("loc/HEK_dpp_enc.bed", fill = T,
                          col.names = c("chr", "start", "end", "cCRE")) %>%
  mutate(cCRE = if_else(cCRE=="", "none", as.character(cCRE))) 

dbts.map$HEK_chh = read.table("loc/HEK_dpp_chh.bed", fill = T,
                          col.names = c("chr", "start", "end", "chromHMM")) %>%
  mutate(chromHMM = if_else(chromHMM=="", "N/A", as.character(chromHMM))) 

dbts.map$HEK_drh = read.table("loc/HEK_dreg_chh.bed", fill = T,
                              col.names = c("chr", "start", "end", "chromHMM")) %>%
  mutate(chromHMM = if_else(chromHMM=="", "N/A", as.character(chromHMM))) 

dbts.map$DPP_enc = read.table("loc/HEK_dpp_enc2.bed", fill = T,
                              col.names = c("chr", "start", "end", "cCRE")) %>%
  mutate(cCRE = if_else(cCRE=="", "none", as.character(cCRE))) 
dbts.map$DRG_enc = read.table("loc/HEK_dreg_enc2.bed", fill = T,
                              col.names = c("chr", "start", "end", "cCRE")) %>%
  mutate(cCRE = if_else(cCRE=="", "none", as.character(cCRE))) 

dbts.hek = bind_rows(
  dbts.map$HEK_chh %>% mutate(method = "deepBTS"),
  dbts.map$HEK_drh %>% mutate(method = "dREG")) %>%
  mutate(chromHMM = factor(chromHMM, levels = c("N/A",
                                        "Ina",
                                        "Rpt",
                                        "Biv",
                                        "Enh",
                                        "Tx",
                                        "TSS")))
dbts.merge = bind_rows(
  dbts.map$HAP %>% select(cCRE) %>% mutate(cell = "HAP1"),
  dbts.map$HEK_enc %>% select(cCRE) %>% mutate(cell = "HEK293"),
) %>%
  mutate(cCRE = factor(cCRE, levels = c("none",
                                        "enhD",
                                        "enhP",
                                        "prom")))

dbts.pr = bind_rows(
  dbts.map$DPP_enc %>% select(cCRE) %>% mutate(method = "deepBTS"),
  dbts.map$DRG_enc %>% select(cCRE) %>% mutate(method = "dREG"),
) %>%
  mutate(cCRE = factor(cCRE, levels = c("none",
                                        "enhD",
                                        "enhP",
                                        "prom")))


library(ggplot2)
# Fraction of cCREs per cell type
pdf("pdf/FigS2D.pdf", width = 2.5, height = 3)
g = ggplot(data = dbts.merge, aes(x = cell)) +
  geom_bar(aes(fill = cCRE), position = "stack", col = "grey60", size = 0.3) +
  scale_fill_manual(values = c("grey", col[c(16,8,5)])) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
print(g)
dev.off()

# cCREs per method
pdf("pdf/Fig2E.pdf", width = 2.5, height = 3)
g = ggplot(data = dbts.pr, aes(x = method)) +
  geom_bar(aes(fill = cCRE), position = "stack", col = "grey60", size = 0.3) +
  scale_fill_manual(values = c("grey", col[c(13,7,4)])) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
print(g)
dev.off()

# Fraction of chromHMM in HEK per method
pdf("pdf/Fig2D.pdf", width = 2.5, height = 3)
g = ggplot(data = dbts.hek, aes(x = method)) +
  geom_bar(aes(fill = chromHMM), position = "stack", col = "grey60", size = 0.3) +
  scale_fill_manual(values = c("grey30", "grey50", col[c(19, 16, 13 ,7,4)])) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
print(g)
dev.off()

# cCRE none/enhD ratio
cCRE.table = dbts.pr %>% group_by(method, cCRE) %>% summarise(count = n())

dbts.pr %>% group_by(method, cCRE) %>% summarise(count = n())
