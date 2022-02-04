setwd("~/Work/labuser/hk572/uPRO/fig/fig7/")
#rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(shades)

# Figs 7A/B - immune relatedgenes
col = rev(colorRampPalette(brewer.pal(11, "RdYlBu")[2:10])(101))
col2 = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(14)
source("rscript/fig7AB.R")

################################
# cCRE analysis (as in Fig 2)
# combined TSS cCRE distribution

col4 = c("grey", colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)[c(13,7,4)])
source("rscript/fig7C.R")


#######################
# TFBS anaylsis (as in Fig 5)
col2 = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)[c(13, 4)]
col3 = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)[c(18, 1)]
source("rscript/fig7D.R")


###########################
# Cell type deconvolution analysis
#
col2 = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(14)
source("rscript/fig7E.R")

###########################################
'Fraction of differentially expressed genes
 after cell type composition normalization'

col2 = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)[c(13, 4)]
source("rscript/fig7E.R")

