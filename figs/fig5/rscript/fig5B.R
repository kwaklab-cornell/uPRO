source("browser/browser.R")
browser.setup(genelist="bed/hg38.refFlat.geneName.bed") # Set up gene list
# Track filename lists
bglist = paste0("bedgraph/",
                c("PBMC_uPRO_1911.pl", "PBMC_uPRO_1911.mn",
                  "PMNL_uPRO_1911.pl", "PMNL_uPRO_1911.mn",
                  "MDM_uPRO_A_2010.pl", "MDM_uPRO_A_2010.mn",
                  "GM18505.pl", "GM18505.mn",
                  "THP1_PRO_U0_1903.pl", "THP1_PRO_U0_1903.mn"),
                ".bedgraph")
#Track colors
color = rep(col2, 5)

# Track heights
heights = rep(0.5, 10)
# Horizontal line positions
hlines = 1:5 * 1
# Track label descriptions
description = c("PBMC", "PMNL", "MDM", "LCL", "THP1")
# Track label positions
label.y = 5:1 * 1 - 0.5
# Setup browser
browser.setup(bedgraph = bglist,
              col = color, heights = heights, hlines = hlines,
              label.y = label.y,
              description = description,
              gene.line = 1)
# Find gene
browser.setpos("chr17", 4300000, 4450000)
#browser.setgene("UBE2G1")
ymax = NULL
ymax = rep(c(700,70,130, 70, 20), each = 2)
browser.read()
ymax = browser.print(filename="pdf/Fig5B.pdf", nbin=400, width=5, height = 3, ymax=ymax)