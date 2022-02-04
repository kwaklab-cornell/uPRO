setwd("fig4/")
#rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(shades)

# Fig4A, read count tables and generate correlation matrix
col = rev(colorRampPalette(brewer.pal(11,"RdYlBu")[2:10])(51))
source("rscript/fig4A.R")


# Fig4B, Differentially expressed TSS, with respect to Cell type specific TSS in Fig3E
col2 = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)[c(13, 4)]
source("rscript/fig4B.R")


# Statistics in the main text
# Differentially expressed in Blood derived cells on average compared to outgroup
diff_frac(ge.ct.bd, "MEAN")
diff_frac(db.ct.bd, "MEAN")
# DE in at least one Blood derived cell compared to outgroup
diff_frac(ge.ct.bd, "UNIT")
diff_frac(db.ct.bd, "UNIT")
# Specific in one cell type within blood derived cells
diff_frac(ge.table, "CTS")
diff_frac(db.table, "CTS")

# Fig4C: heatmap of DE TSSs
col = rev(colorRampPalette(brewer.pal(11,"RdYlBu")[2:10])(101))
col2 = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)[c(13, 4)]
source("rscript/fig4C.R")

# Fig 4D: PCA analysis
col = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)
source("rscript/fig4D.R")

# FigS4
col = rev(colorRampPalette(brewer.pal(11,"RdYlBu")[2:10])(51))
col2 = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)[c(16, 5)]
col3 = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)
source("rscript/figS4.R")



  
