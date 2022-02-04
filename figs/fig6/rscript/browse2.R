source("browser/browser.R")
library(RColorBrewer)

# Set up gene list
browser.setup(genelist="geneAnnotations/gencode.v26.uniq.bed") # Set up gene list

# Track filename lists
bglist = paste0("bedgraph/",
		c("pbmc.pl", "pbmc.mn",
		  "pmnl.pl", "pmnl.mn"),
		".bedgraph")

# Track colors
color.pm = colorRampPalette(brewer.pal(11, "Spectral"))(20)
color = c(color.pm[c(5,19)],
	  color.pm[c(5,19)])

# Track heights
heights = c(rep(0.75, 4))

# Horizontal line positions
hlines = c(3, 1.5)

# Track label descriptions
description = c("PBMC", "PMNL")

# Track label positions
label.y = c(2.25, 0.75)

# Setup browser
browser.setup(bedgraph = bglist,
              col = color, heights = heights, hlines = hlines,
              label.y = label.y,
              description = description)

tc = c(
8025410 + 7810593,
971141 + 825799 ) / 1000000
ymax = rep(tc, each = 2) * c(20)

# Find gene
browser.setpos("chr19", 43500000, 45500000)
#ymax = NULL
browser.read()
ymax = browser.print(filename="pdf/Fig5C.pdf", nbin=400, width=8, height = 2, ymax=ymax)
