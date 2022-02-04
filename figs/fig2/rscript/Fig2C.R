setwd("~/Work/labuser/hk572/uPRO/fig/fig2/")

# Browser view of various transcription data
source("browser/browser.R")
browser.new()

# Set up gene list
browser.setup(genelist="bed/hg38.refFlat.geneName.bed")

# Track filename lists
bglist = c("bg/HEK.pl.bedgraph", "bg/HEK.mn.bedgraph",
           "bg/nuPRO.pl.bedgraph", "bg/nuPRO.mn.bedgraph",
           "bed/HAP1.dREG.bedgraph",
           "bed/HAP1_uPRO_1910.dBTS.bedgraph")

# Track colors
library(RColorBrewer)
color = c(col[c(4,17)],
          col[c(2,19)],
          "grey50",
          "grey30")

# Track heights
heights = c(1, 1,
            1, 1,
            0.75,
            0.75)

# Horizontal line positions
hlines = c(0.75,
           1.5,
           3.5,
           5.5)

# Track label descriptions
description = c("HEK293",
                "HAP1",
                "dREG",
                "deepBTS")

# Track label positions
label.y = c(4.5,
            2.5,
            1.125,
            0.375)

# Setup browser
browser.setup(bedgraph = bglist,
              col = color, heights = heights, hlines = hlines,
              label.y = label.y,
              description = description,
              gene.line = 1)

# Find gene
if(F) {
browser.setgene("HIVEP2", mar = 0.4)
ymax = c(18000, 18000,
       8000, 8000,
       50000, 500,
       7000, 4000, 5000)
ymax = NULL
browser.read()
browser.print(filename="Fig2C.pdf", nbin=400,
              width = 8, height = 3, ymax=ymax)
}

# Find position
browser.setpos("chr2",180143350-250000, 180148150+250000)
browser.read()
ymax = c(250, 250, 250, 250, 35, 35)
browser.print(filename="pdf/Fig2C.pdf" , nbin=400,
              width = 8, height = 3, ymax=ymax)

