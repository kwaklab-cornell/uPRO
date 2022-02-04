
aurocfiles = c("roc/auroc.GM19238.dreg.txt",
               "roc/GM19238.auroc.0607-1-1.txt")
labs = c("dREG", "deepBTS")

read_roc = function() {
	data = NULL
	for(i in 1:length(aurocfiles)) {
	  data = rbind(data,
		       data.frame(Method = labs[i],
				  read.table(aurocfiles[i], header = T)))
	  data = rbind(data, data.frame(Method = labs[i],
				  cutoff = c(-Inf, Inf),
				  pos_evl = c(0, 1),
				  tot_evl = c(1, 1),
				  pos_std = c(1, 0),
				  tot_std = c(1, 1)))
	}
	data = data %>%
	  mutate(FPR = 1 - pos_evl/tot_evl, TPR = pos_std/tot_std)
	return(data)
}

data = read_roc()

calc_auc = function(data) {
  d = data.frame(x = data$FPR, y = data$TPR) %>%
    na.omit %>%
    arrange(x) %>%
    mutate(w = c(x[-1],1) - x)
  return(sum(d$y * d$w))
}

aucs = data.frame(x = rep(0.5, 2),
                  y = c(0.175, 0.1),
                  label = c("dREG", "deepBTS"))

aucs= aucs %>%
  mutate(label = paste0(label, ": ", 
                        sapply(label ,function(x) round(calc_auc(data %>% filter(Method == x)), digits = 3))))

pdf("pdf/Fig2B.pdf", width = 4, height = 3)
g = ggplot(data, aes(x = FPR, y = TPR, col = Method)) +
  geom_vline(xintercept = 0.1, linetype = "dotted", size = 0.5) +
  geom_step(size = 0.8) +
  theme_bw() +
  scale_color_manual(values = col) +
  geom_text(data = aucs, aes(x = x, y = y, hjust = 0, vjust = 1,
                label = label), col = "black")
print(g)
dev.off()


plot_roc = function(roc_file, label) {
  aurocfiles[3] <<- roc_file
  labs[3] <<- label
  data <<- read_roc()
  ggplot(data, aes(x = FPR, y = TPR, col = Method)) +
    geom_vline(xintercept = 0.1, linetype = "dotted", size = 0.5) +
    geom_step()
  print(paste("dREG =" ,calc_auc(data %>% filter(Method == "dREG"))))
  print(paste(label, "=" ,calc_auc(data %>% filter(Method != "dPRO07-1" & Method != "dREG"))))
}

