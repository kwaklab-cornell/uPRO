# function to remove overlapping spots
reduce.scatter = function(data, x, y, grid = 100)
{
  xlim = range(x)
  ylim = range(y)
  den.mat = matrix(0, ncol = grid + 1, nrow = grid + 1)
  x.grid = floor( ( x - xlim[1] ) / ( xlim[2] - xlim[1] ) * grid + 1)
  y.grid = floor( ( y - ylim[1] ) / ( ylim[2] - ylim[1] ) * grid + 1)
  data2 = data %>% mutate(incl = F)
  n = nrow(data2)
  for(i in 1:n) {
    if( den.mat[ x.grid[i], y.grid[i] ] == 0 ) {
      data2$incl[i] = T
      den.mat[ x.grid[i], y.grid[i] ] = 1
    }
  }
  return(data2 %>% filter(incl == T) %>% select(-incl))
}

# Correlation plot function
plot.cor = function(data, label_x = '', label_y = '',
                    red.level = 400, base = 10,
                    xlim = NULL, ylim = NULL,
                    diag = TRUE, cor = FALSE)
{
  minion.colors = c("#8064A260","#1F497D60","#77933C60","#F0C01060",
                    "#D0701060","#95373560") 
  data2 = data %>%
    filter(x > 0 & y > 0)
  plot.rr = data2 %>%
    mutate(density = densCols(log(x)/log(base), log(y)/log(base),
                              nbin=512, bandwidth = 0.1,
                              colramp = colorRampPalette(minion.colors)))
  if(red.level > 0) plot.rrr = reduce.scatter(plot.rr,
                            log(plot.rr$x)/log(base),
                            log(plot.rr$y)/log(base),
                            red.level)
  else plot.rrr = plot.rr
  fitdata = data.frame(x = log(data2$x)/log(base), y = log(data2$y)/log(base))
  pc1 = prcomp(fitdata)$rotation
  slope = pc1[2,1]/pc1[1,1]
  intercept = mean(fitdata$y) - slope * mean(fitdata$x)
  corr = round(cor(fitdata$x, fitdata$y), digits = 3)
  plot.obj = ggplot(plot.rrr, aes(y = log(y)/log(base), x = log(x)/log(base))) + 
    geom_point(aes(col = density), size=1) + scale_color_identity() +
    theme_bw() + 
    labs(x = label_x, y = label_y)
  if(!is.null(xlim)) plot.obj = plot.obj + xlim(xlim)
  if(!is.null(ylim)) plot.obj = plot.obj + ylim(ylim)
  if(diag) plot.obj = plot.obj + geom_abline(slope = slope, intercept = intercept)
  if(cor) plot.obj = plot.obj +
	  geom_text(aes(x = Inf, y = -Inf, hjust = 1.1, vjust = -0.25,
			label = corr, fontface = "italic"))
  return(plot.obj)
}

# Scatter plot function
plot.scatter = function(data, label_x = '', label_y = '',
                    red.level = 400, fit = "none",
                    xlim = NULL, ylim = NULL)
{
  minion.colors = c("#8064A260","#1F497D60","#77933C60","#F0C01060",
                    "#D0701060","#95373560") 
  plot.rr = data %>%
    mutate(density = densCols(x, y,
                              nbin=1024, bandwidth = 0.1,
                              colramp = colorRampPalette(minion.colors)))
  if(red.level > 0) plot.rrr = reduce.scatter(plot.rr,
                                              plot.rr$x,
                                              plot.rr$y,
                                              red.level)
  else plot.rrr = plot.rr
  if(fit == "lowess") {
    ls1 = lowess(y ~ x, data)
  }
  
  plot.obj = ggplot(plot.rrr, aes(y = y, x = x)) + 
    geom_point(aes(col = density), size=1) + scale_color_identity() +
    theme_bw() + 
    labs(x = label_x, y = label_y)
  if(!is.null(xlim)) plot.obj = plot.obj + xlim(xlim)
  if(!is.null(ylim)) plot.obj = plot.obj + ylim(ylim)
  if(fit == "lowess") plot.obj = plot.obj
  return(plot.obj)
}

