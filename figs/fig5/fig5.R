setwd("fig5/")
#rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(shades)

# Fig 5A: Load and process correlation matrices, generate network star shape
col = rev(colorRampPalette(brewer.pal(11,"RdYlBu")[2:10])(100))
col2 = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)[c(3, 17)]
source("rscript/fig5A.R")

# Fig 5B example genome browser at CYP11A1 gene
col2 = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)[c(3, 18)]
source("rscript/fig5B.R")

# Fig 5C correlation scatterplots
col4 = brightness(colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)[c(1,4,7,16)], 0.85)
source("rscript/fig5C.R")

# Fig5D : properties of TF-enhancer coexpression - which TFs are better coexpressed?
col2 = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)[c(13, 4)]
col3 = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)[c(18, 1)]
source("rscript/fig5D.R")

# Fig5E : properties of enhancer-TF coexpression - which enhancers are better co-expressed?
col2 = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)[c(13, 4)]
col3 = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)[c(18, 1)]
source("rscript/fig5E.R")

# Fig5F : properties of enhancer-gene coexpression - distance effect
source("rscript/fig5F.R")
