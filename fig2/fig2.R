# Fig 2 
setwd("fig2/")
#rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(shades)

# Fig 2B
col = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)[c(16, 5)]
source("rscript/Fig2B.R")

# Fig 2C
col = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)
source("rscript/Fig2C.R")

# Fig 2DE
# run sh/fig2D.sh first

col = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)
source("rscript/Fig2DE.R")

# Other sunoke counts

# Cell type specific enhancers example
rc.table = read.table("readcounts.txt", header = T)
rc.enh = rc.table %>%
  filter(grepl("enh", id)) %>%
  arrange(-HAP1_PRO_1904)
  
# Fractions
dbts.hek %>%
  filter(method == "dREG") %>% 
  filter(chromHMM == "N/A" | chromHMM == "Ina" | chromHMM == "Rpt") %>%
  nrow /
  dbts.hek %>% filter(method == "dREG") %>% nrow 
