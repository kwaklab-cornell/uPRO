################################
# cCRE analysis (Fig 2)
# combined TSS cCRE distribution
ccre = 
  read.table("bed/WB_cCRE.bed", fill = T,
             col.names = c("chr", "start", "end", "cCRE")) %>%
  mutate(cCRE = if_else(cCRE=="", "none", as.character(cCRE))) %>%
  mutate(sample = "Whole blood") %>%
  mutate(cCRE = factor(cCRE,
                       levels = c("none",
                                  "enhD",
                                  "enhP",
                                  "prom")))


# ccre count plot
pdf("pdf/Fig7C.pdf", width = 2, height = 2.5)
g = ggplot(data = ccre, aes(x = sample)) +
  geom_bar(aes(fill = cCRE), position = "stack", col = "grey60") +
  scale_fill_manual(values = col4) +
  theme_bw()
print(g)
dev.off()

# ccre percents
ccre.table = ccre %>%
  group_by(cCRE) %>%
  summarise(count = n()) %>%
  mutate(percent = count / sum(count) * 100)