###########################
# Cell type deconvolution analysis
#

library(dplyr)
library(tidyr)

tr.s = read.table("deconv/sigExp.txt", header = T, stringsAsFactors = F)
pchro = ge.rc %>%
  separate(id, c("pos", "name"), sep = ";") %>%
  select(-pos) %>%
  distinct(name, .keep_all = TRUE)

pchro.name = pchro

pchro.all = pchro.name %>%
  inner_join(tr.s, by = "name")

pchro.plot = pchro.all
pchro.plot[,2:13] = pchro.plot[,2:13]/pchro.plot[,2]
pchro.plot = pchro.plot %>%
  select(-P1, -count, -r)

######################################################
pchro.avsig = pchro.plot %>%
  select(-PBMC, -PMNL)%>%
  gather(sample, val, -name, -sig) %>%
  filter(is.finite(val) & val > 0) %>%
  group_by(sig, sample) %>%
  summarise(mean = 2^mean(log2(val)))%>%
  ungroup() %>%
  spread(sig, mean) %>%
  mutate(PBMC = PBMC / NONE, PMNL = PMNL / NONE) %>%
  select(-NONE) %>%
  mutate(sample = paste0("mean_", sample)) %>%
  gather(sig, mean, -sample) %>%
  spread(sample, mean) 

pchro.sesig = pchro.plot %>%
  select(-PBMC, -PMNL) %>%
  gather(sample, val, -name, -sig) %>%
  filter(is.finite(val) & val > 0) %>%
  group_by(sig, sample) %>%
  summarise(se = 2^sd(log2(val)/sqrt(n()-1))) %>%
  ungroup() %>%
  spread(sig, se) %>%
  select(-NONE) %>%
  mutate(sample = paste0("se_", sample)) %>%
  gather(sig, se, -sample) %>%
  spread(sample, se) 

pchro.sig = pchro.plot %>%
  filter(sig != "NONE") %>%
  arrange(sig)


##################################################

refFrac = data.frame(sig = c("PBMC", "PMNL"),
                     rfr = c(0.35, 0.65))

dcp.plot = pchro.avsig %>%
  gather(sample, mean, -sig) %>%
  inner_join(refFrac, by = "sig") %>%
  mutate(fraction =  mean * rfr) %>%
  select(sig, sample, fraction) %>%
  spread(sample, fraction) %>%
  mutate_if(is.numeric, funs(./sum(.))) %>%
  gather(sample, fraction, -sig) %>%
  separate(sample, c("type", "indiv"), extra = "merge") %>%
  select(sig, indiv, fraction) %>%
  inner_join(pchro.sesig %>%
               gather(sample, se, -sig) %>%
               separate(sample, c("type", "indiv"), extra = "merge") %>%
               select(sig, indiv, se),
             by = c("sig", "indiv")) %>%
  mutate(fr_ue = se * fraction,
         fr_le = fraction / se) %>%
  mutate(indiv = factor(indiv,
                        levels = c("P1r", "P1n","P1nh",
                                   "P2", "P6", "P33",
                                   "P45", "P52", "P53")))

# Plot barplot with error bars
pdf("pdf/Fig7E.pdf", width = 4, height =3)
par(mar = c(2, 2, 1, 1), mgp = c(1.5, 0.5, 0))
g = ggplot(data = dcp.plot, aes(x = sig, y= fraction, fill = indiv)) +
  geom_bar(stat = "identity",
           position = position_dodge(0.9),
           col = "grey60",
           size = 0.25) +
  geom_errorbar(aes(ymin = fr_le, ymax = fr_ue),
                width = 0.4, position=position_dodge(0.9),
                col = "grey60", size = 0.3) +
  scale_fill_manual(values = col2[c(2:7, 9, 11:12)]) +
  ylab("Cell type fraction") +
  xlab("Cell type") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(legend.key = element_rect())
print(g)
dev.off()
