# TF-enhancer-gene coexpression network
setwd("tegcn/")
source("scripts/tegcn.R")

# load and calculation coexpression coefficient in
# db.gn and tf.en data
tcn = make_tegcn()

# Plot network on chr 17
tcn.chr17.net = make_2tfnet(tcn, tf1 = "E47", tf2 = "USF", chr = "chr17")
plot_tfnet(tcn.chr17.net)

# Tf-enhancer correlations
tf.cor = tegcn_tfcor(tegcn = tcn)
plot_cognate(tf.cor)

# TFs that are significantly more correlated to cognate than non-cognate
tf.sig = tfcor_sig(tf.cor)
plot_guitar(tf.cor, n = 20)

# Distance plot
plot_cordist(tcn)


