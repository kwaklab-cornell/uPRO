# combined TSS cCRE distribution
read.ccre = function(name) {
  return(
    read.table(paste0("ccre/", name, ".bed"), fill = T,
               col.names = c("chr", "start", "end", "cCRE")) %>%
      mutate(cCRE = if_else(cCRE=="", "none", as.character(cCRE))) %>%
      mutate(sample = name)) %>%
    mutate(cCRE = factor(cCRE,
                         levels = c("none",
                                    "enhD",
                                    "enhP",
                                    "prom")))
}

combined.ccre = read.ccre("combined")
combined.table = combined.ccre %>%
  group_by(cCRE) %>%
  summarise(count = n()) %>%
  mutate(percent = count / sum(count) * 100)

# Every TSS cCRE distribution
samples = c(
  "HAP1_PRO_1904", "HAP1_uPRO_1910",
  "HEK_PRO_DMSO1_1701", "HEK_PRO_DMSO2_1701", "HEK_PRO_THAP1_1701", "HEK_PRO_THAP2_1701",
  "HeLa_PRO_1800", "HeLa_PRO_WN_1903", "HeLa_uPRO_ARF_2003", "HeLa_uPRO_SSK_2003",
  "LCL_GM18505_1500", "LCL_GM18517_1500", "LCL_GM18517r_1500", "LCL_GM18520_1500", "LCL_GM18520r_1500",
  "LCL_GM18522_1500", "LCL_GM18522r_1500", "LCL_GM19099r_1500", "LCL_GM19193_1500", "LCL_GM19193r_1500",
  "LCL_GM19222_1500", "LCL_GM19222r_1500", "LCL_GM19238_1500", "LCL_GM19238r_1500",
  "LCL_GM19239_1500", "LCL_GM19239r_1500",
  "MDM_uPRO_A_2010", "MDM_uPRO_B_2010", "MDM_uPRO_C_2010",
  "PBMC_uPRO_1011", "PMNL_uPRO_1911",
  "THP1_PRO_D0_1903", "THP1_PRO_D1_1903", "THP1_PRO_D2_1903", "THP1_PRO_D4_1903",
  "THP1_PRO_U0_1903", "THP1_PRO_U2_1903", "THP1_PRO_U4_1903"
)

ccre.table = data.frame()
for(name in samples) {
  ccre.table = ccre.table %>%
    bind_rows(read.ccre(name))
}
ccre.table = ccre.table %>%
  mutate(cCRE = factor(cCRE,
                       levels = c("none",
                                  "enhD",
                                  "enhP",
                                  "prom")))
ccre.table = ccre.table %>%
  mutate(sample = gsub("_[0-9]+$", "", sample))

ccre.split = ccre.table %>%
  group_by(sample, cCRE) %>%
  summarise(count = n()) %>%
  mutate(percent = count / sum(count) * 100) %>%
  mutate(total = sum(count))

# Fig 3A, precentage of each cCRE elements
pdf("pdf/Fig3A.pdf", width = 2, height = 3)
g = ggplot(data = ccre.split) +
  geom_boxplot(aes(x = cCRE, y = percent, fill = cCRE), col = "grey60") +
  scale_fill_manual(values = c("grey", col[c(13,7,4)])) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_discrete(limits = c("prom", "enhP", "enhD", "none")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
print(g)
dev.off()

# Distribution of TSS numbers
ccre.split %>%
  group_by(cCRE) %>%
  summarize(mean = mean(percent))

# TSS discovery per read count
total.rc = read.table("readcounts/totalreadcounts.txt", col.names = c("sample", "readcount")) %>%
  bind_cols(ccre.split %>%
              filter(cCRE == "none"))

# Fig 3B, read count vs TSS
pdf("pdf/Fig3B.pdf", width = 3, height = 3)
g = ggplot(data = total.rc, aes(x = readcount/1000000, y = total)) +
  geom_point() +
  theme_bw() +
  xlab("Readcount (million)") +
  ylab("deepBTS TSSs") +
  ylim(0, 25000)
print(g)
dev.off()
