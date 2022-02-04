source("browser/browser.R")
library(RColorBrewer)

# Set up gene list
browser.setup(genelist="geneAnnotations/gencode.v26.name.bed") # Set up gene list

# Track filename lists
bglist = paste0("HAP1/",
		c("nuPRO.pl", "nuPRO.mn",
		  "PRO.pl", "PRO.mn",
		  "HEK.pl", "HEK.mn"),
		".bedgraph")

# Track colors
color = c(col[c(2,19)],
	  col[c(2,19)],
	  col[c(4,17)])

# Track heights
heights = c(rep(0.75, 6))

# Horizontal line positions
hlines = c(4.5, 3, 1.5)

# Track label descriptions
description = c("HAP1\nuPRO",
                "HAP1\nPRO-seq",
		"HEK293")

# Track label positions
label.y = c(3.75, 2.25, 0.75)

# Setup browser
browser.setup(bedgraph = bglist,
              col = color, heights = heights, hlines = hlines,
              label.y = label.y,
              description = description)

# Find gene
browser.setpos("chr1",47200000,47325000)
#browser.setgene("TAL1")
ymax = c(50, 50, 50, 50, 100, 100)
browser.read()
ymax = browser.print(filename="pdf/Fig1A.pdf", nbin=400, width=6, height = 3, ymax=ymax)



