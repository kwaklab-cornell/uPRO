setwd("~/Work2/hk572/tiQTL_browser/")    # Change this to your work directory location
source("browser.R")
source("viewQTL.R")

show.procap("rs17149633")
tiQTL.list = read.table("data/list/Table.9a.high.causal.tiQTL.txt", header = 1)
