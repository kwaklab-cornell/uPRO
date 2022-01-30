source("browser/browser.R")
library(RColorBrewer)

# Set up gene list
browser.setup(genelist="geneAnnotations/gencode.v26.name.bed") # Set up gene list

# Track filename lists
bglist = paste0("bedgraph/",
		c("hk.pl", "hk.mn",
		  "sl.pl", "sl.mn",
		  "pbmc.pl", "pbmc.mn",
		  "pmnl.pl", "pmnl.mn"),
		".bedgraph")

# Track colors
color.pm = colorRampPalette(brewer.pal(11, "Spectral"))(20)
color = c(color.pm[c(3,18)],
	  color.pm[c(3,18)],
	  color.pm[c(5,19)],
	  color.pm[c(5,19)])

# Track heights
heights = c(rep(0.75, 8))

# Horizontal line positions
hlines = c(6, 4.5, 3, 1.5)

# Track label descriptions
description = c("pChRO\nIndiv 1",
                "pChRO\nIndiv 2",
		"PBMC", "PMNL")

# Track label positions
label.y = c(5.25, 3.75, 2.25, 0.75)

# Setup browser
browser.setup(bedgraph = bglist,
              col = color, heights = heights, hlines = hlines,
              label.y = label.y,
              description = description)

tc = c(
3753359 + 3433510,
1058066 + 906826,
8025410 + 7810593,
971141 + 825799 ) / 1000000
ymax = rep(tc, each = 2) * c(20, 10)

# Find gene
browser.setpos("chr1", 21499917, 21800345)
# ymax = NULL
#ymax = c(250, 250, 70, 70, 40, 40)
browser.read()
ymax = browser.print(filename="pdf/Fig5A.pdf", nbin=200, width=6, height = 3.5, ymax=ymax)
