# Hierarchical TF-enhancer-gene coexpression network (TEGCN) analysis for uPRO

```
setwd("tegcn/")
source("scripts/tegcn.R")
```

load and calculation coexpression coefficient in
```
tcn = make_tegcn()
```

Plot network on chr 17
```
tcn.chr17.net = make_2tfnet(tcn, tf1 = "E47", tf2 = "USF", chr = "chr17")
plot_tfnet(tcn.chr17.net)
```
![plot](./tegcn/plot/chr17_tfnet.png)


TFs-enhancer correlations that are significantly more correlated to cognate than non-cognate
```
tf.cor = tegcn_tfcor(tegcn = tcn)
tf.sig = tfcor_sig(tf.cor)
plot_guitar(tf.cor, n = 20)
```

Distance plot
```
plot_cordist(tcn)
```

