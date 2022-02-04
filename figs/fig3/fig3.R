setwd("fig3/")
#rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(shades)

# Fig 3AB
col = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)
source("rscript/fig3AB.R")


# Fig 3C, after running cellstypespecific.sh, cell type specific enhancer counts
col = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)
col2 = col[c(4, 16)]
source("rscript/fig3CEDF.R")

# Fig S3
col = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)
source("rscript/figS3.R")