# Fig 1 
setwd("fig1/")
#rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(shades)

# Fig 1A
col = colorRampPalette(brewer.pal(11,"RdYlBu")[1:11])(20)
source("rscript/Fig1A.R")

# Fig 1BCDE
col = rev(colorRampPalette(brewer.pal(11,"RdYlBu")[2:10])(51))
source("rscript/Fig1BCDE.R")
