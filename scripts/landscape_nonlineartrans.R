rm(list = ls())
set.seed(0)

# load auxiliary functions
source('./ecoFunctions.R')
library(tidyr)
library(CEGO)
library(ggExtra)

# tunable parameters
n_mut <- 6
mycolors <- c('#939598', '#be1e2d', '#85c441', '#d68f28')

# make genotype names
genots <- lapply(0:n_mut, FUN = function(i) t(combn(n_mut, i)))
genots <- lapply(genots, FUN = function(x) sapply(1:nrow(x),
                                                  FUN = function(i) paste(x[i, ], collapse = ',')))
genots <- unlist(genots)

# make edges of fitness graph
nMut <- function(genot) sapply(genot,
                               FUN = function(genot_i) length(strsplit(genot_i, split = ',')[[1]]))

isDescendant <- function(this_genot, of_this_genot) {
  
  this_genot <- strsplit(this_genot, split = ',')[[1]]
  of_this_genot <- strsplit(of_this_genot, split = ',')[[1]]
  
  return(all(of_this_genot %in% this_genot))
  
}

edges <- data.frame(source = character(0),
                    target = character(0),
                    source.nmut = numeric(0),
                    target.nmut = numeric(0))

for(s in genots) {
  
  t <- genots[sapply(genots,
                     isDescendant,
                     of_this_genot = s) & nMut(genots) == nMut(s)+1]
  if(length(t)) {
    edges <- rbind(edges,
                   data.frame(source = s,
                              target = as.character(t),
                              source.nmut = as.numeric(nMut(s)),
                              target.nmut = as.numeric(nMut(s)) + 1))
  }
  
}

edges <- cbind(edge_id = paste('edge_', 1:nrow(edges), sep = ''),
               edges)

# plot landscape (wrapper function)
plotGraph <- function(landscape, save.plot = F) {

  df <- cbind(edges,
              source.f = setNames(landscape$f, landscape$genot)[edges$source],
              target.f = setNames(landscape$f, landscape$genot)[edges$target])
  df$source.f[is.na(df$source.f)] <- landscape$f[landscape$genot == '']
  
  if ('color' %in% colnames(landscape)) {
    df <- merge(df, landscape[, c('genot', 'color')], by.x = 'target', by.y = 'genot')
  } else {
    df$color <- 'A'
  }
  df <- df[, c('edge_id', 'source', 'target', 'source.nmut', 'target.nmut', 'source.f', 'target.f', 'color')]
  
  dfx <- gather(df[, c(1, 4, 5)], position, nmut, source.nmut:target.nmut)
  dfx$position <- setNames(c('source', 'target'), c('source.nmut', 'target.nmut'))[dfx$position]
  
  dfy <- gather(df[, c(1, 6, 7)], position, f, source.f:target.f)
  dfy$position <- setNames(c('source', 'target'), c('source.f', 'target.f'))[dfy$position]
  
  dfxy <- merge(dfx, dfy, by = c('edge_id', 'position'))
  
  df <- merge(dfxy, df[, c('edge_id', 'color')], by = 'edge_id')
  
  myplot <-
    ggplot(df, aes(x = nmut, y = f, group = edge_id, color = color)) +
      geom_abline(slope = 0,
                  intercept = landscape$f[landscape$genot == ''],
                  color = '#d1d3d4') +
      geom_line() +
      scale_x_continuous(name = '# mutations',
                         breaks = 0:n_mut,
                         labels = as.character(0:n_mut)) +
      scale_y_continuous(name = expression(italic(F)),
                         breaks = pretty_breaks(n = 3),
                         expand = c(0.05, 0.05)) +
      scale_color_manual(values = setNames(mycolors, LETTERS[1:length(mycolors)])) +
      theme_bw() +
      theme(aspect.ratio = 0.6,
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = 'none',
            axis.title = element_text(size = 18),
            axis.text = element_blank(), #axis.text = element_text(size = 16),
            axis.ticks.length = unit(0, 'cm')) + #axis.ticks.length = unit(-0.15, 'cm')) +
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)
  
  if (save.plot != F) {
    ggsave(myplot,
           file = paste('../plots/', save.plot, '.pdf', sep = ''),
           dpi = 600,
           width = 60,
           height = 50,
           units = 'mm')
  }
  
  return(myplot)

}

# plot FEEs (wrapper function)
plotFEEs <- function(landscape, save.plot = F) {
  
  gedf <- makeGEdata(landscape[, c(1, 2)])
  dfgenots <- sapply(1:nrow(gedf),
                     FUN = function(i) paste(sort(c(strsplit(gedf$background[i], split = ',')[[1]], gedf$knock_in[i])), collapse = ','))
  if('color' %in% colnames(landscape)) {
    gedf$color <- setNames(landscape$color, landscape$genot)[dfgenots]
  } else {
    gedf$color <- 'A'
  }
  
  gedf$knock_in <- paste('mut.', gedf$knock_in)
  
  myplot <-
    ggplot(gedf, aes(x = background_f, y = d_f, group = knock_in, color = color)) +
      geom_abline(slope = 0,
                  intercept = 0,
                  color = '#d1d3d4') +
      geom_point(cex = 2) +
      geom_smooth(method = 'lm',
                  formula = y~x,
                  se = F,
                  fullrange = T,
                  color = 'black') +
      facet_wrap(~knock_in) +
      scale_x_continuous(name = expression(paste(italic(F), ' (background)', sep = '')),
                         breaks = pretty_breaks(n = 3),
                         expand = c(0.05, 0.05)) +
      scale_y_continuous(name = expression(paste(Delta, italic(F), sep = '')),
                         breaks = pretty_breaks(n = 2),
                         expand = c(0.05, 0.05)) +
      scale_color_manual(values = setNames(mycolors, LETTERS[1:length(mycolors)])) +
      theme_bw() +
      theme(aspect.ratio = 0.6,
            panel.grid = element_blank(),
            panel.border = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(size = 16,
                                      hjust = 0,
                                      face = 'italic'),
            legend.position = 'none',
            axis.title = element_text(size = 18),
            axis.title.x = element_text(vjust = 0),
            axis.text = element_blank(), #axis.text = element_text(size = 16),
            axis.ticks.length = unit(0, 'cm')) + #axis.ticks.length = unit(-0.15, 'cm')) +
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)
  
  if (save.plot != F) {
    ggsave(myplot,
           file = paste('../plots/', save.plot, '.pdf', sep = ''),
           dpi = 600,
           width = 80,
           height = 65,
           units = 'mm')
  }
  
  return(myplot)
  
}

# plot non-linear transformation (wrapper function)
plotNLT <- function(landscape, nlt = function(x) x, save.plot = F) {
  
  myplot <-
    ggplot(data.frame(x = (landscape$f - min(landscape$f))/(max(landscape$f) - min(landscape$f)),
                      y = y),
           aes(x = x, y = y)) +
    geom_point(alpha = 0.25,
               cex = 2,
               shape = 16) +
    geom_line(data = data.frame(x = seq(0, 1, length.out = 200),
                                y = nlt(seq(0, 1, length.out = 200))),
              aes(x = x, y = y)) +
    scale_x_continuous(name = expression(lambda)) +
    scale_y_continuous(name = expression(paste(italic(f), ' (', lambda, ')', sep = ''))) +
    theme_bw() +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 16,
                                    hjust = 0,
                                    face = 'italic'),
          legend.position = 'none',
          axis.title = element_text(size = 18),
          axis.title.x = element_text(vjust = 0),
          axis.text = element_blank(), #axis.text = element_text(size = 16),
          axis.ticks.length = unit(0, 'cm')) #axis.ticks.length = unit(-0.15, 'cm'))
  
  myplot <- ggMarginal(myplot, size = 8)
  print(myplot)
  
  if (save.plot != F) {
    ggsave(myplot,
           file = paste('../plots/', save.plot, '.pdf', sep = ''),
           dpi = 600,
           width = 40,
           height = 40,
           units = 'mm')
  }
  
}

### LANDSCAPES

# 1: additive landscape
s <- setNames(sample(c(rnorm(3, mean = 0.1, sd = 0.05),
                       rnorm(3, mean = -0.1, sd = 0.05))),
              1:n_mut)
landscape <- data.frame(genot = genots,
                        f = sapply(genots,
                                   FUN = function(genot_i) 1 + sum(s[strsplit(genot_i, split = ',')[[1]]])))

plotGraph(landscape, save.plot = 'test1')
plotFEEs(landscape, save.plot = 'test2')

# save additive landscape
landscape0 <- landscape

# 2: concave non-linear transformation
nlt <- function(x) (1-exp(-5*x))/(1-exp(-5)) # convention: min. for x=0, max. for x=1

landscape <- landscape0

y <- nlt((landscape$f - min(landscape$f))/(max(landscape$f) - min(landscape$f)))
landscape$f <- y*(max(landscape$f) - min(landscape$f)) + min(landscape$f)

plotNLT(landscape0, nlt, save.plot = 'NLT_concave_transformation')
plotGraph(landscape, save.plot = 'NLT_concave_graph')
plotFEEs(landscape, save.plot = 'NLT_concave_fees')

# 3: convex non-linear transformation
nlt <- function(x) (exp(5*x)-1)/(exp(5)-1) # convention: min. for x=0, max. for x=1

landscape <- landscape0

y <- nlt((landscape$f - min(landscape$f))/(max(landscape$f) - min(landscape$f)))
landscape$f <- y*(max(landscape$f) - min(landscape$f)) + min(landscape$f)

plotNLT(landscape0, nlt, save.plot = 'NLT_convex_transformation')
plotGraph(landscape, save.plot = 'NLT_convex_graph')
plotFEEs(landscape, save.plot = 'NLT_convex_fees')

# 4: concave non-monotonic transformation
nlt <- function(x) 1 - 4*(x - 0.5)^2 # convention: min. for x=0, max. for x=1

landscape <- landscape0

y <- nlt((landscape$f - min(landscape$f))/(max(landscape$f) - min(landscape$f)))
landscape$f <- y*(max(landscape$f) - min(landscape$f)) + min(landscape$f)

plotNLT(landscape0, nlt, save.plot = 'NLT_concaveNM_transformation')
plotGraph(landscape, save.plot = 'NLT_concaveNM_graph')
plotFEEs(landscape, save.plot = 'NLT_concaveNM_fees')

# 5: convex non-monotonic transformation
nlt <- function(x) 4*(x - 0.5)^2 # convention: min. for x=0, max. for x=1

landscape <- landscape0

y <- nlt((landscape$f - min(landscape$f))/(max(landscape$f) - min(landscape$f)))
landscape$f <- y*(max(landscape$f) - min(landscape$f)) + min(landscape$f)

plotNLT(landscape0, nlt, save.plot = 'NLT_convexNM_transformation')
plotGraph(landscape, save.plot = 'NLT_convexNM_graph')
plotFEEs(landscape, save.plot = 'NLT_convexNM_fees')

# 6: convex to concave transformation
nlt <- function(x) {
  
  k <- 10
  
  L <- (1 + cosh(k/2))/sinh(k/2)
  C <- L/(1 + exp(k/2))
  
  return(L/(1 + exp(-k*(x-0.5))) - C)
  
   
}

landscape <- landscape0

y <- nlt((landscape$f - min(landscape$f))/(max(landscape$f) - min(landscape$f)))
landscape$f <- y*(max(landscape$f) - min(landscape$f)) + min(landscape$f)

plotNLT(landscape0, nlt, save.plot = 'NLT_logistic-convex-to-concave_transformation')
plotGraph(landscape, save.plot = 'NLT_logistic-convex-to-concave_graph')
plotFEEs(landscape, save.plot = 'NLT_logistic-convex-to-concave_fees')

# 7: concave to convex transformation
nlt <- function(x) {
  
  k <- 10
  
  L <- (1 + cosh(k/2))/sinh(k/2)
  C <- L/(1 + exp(k/2))
  
  return(0.5 - (1/k)*log(L/(x+C) - 1))
  
  
}

landscape <- landscape0

y <- nlt((landscape$f - min(landscape$f))/(max(landscape$f) - min(landscape$f)))
landscape$f <- y*(max(landscape$f) - min(landscape$f)) + min(landscape$f)

plotNLT(landscape0, nlt, save.plot = 'NLT_logistic-concave-to-convex_transformation')
plotGraph(landscape, save.plot = 'NLT_logistic-concave-to-convex_graph')
plotFEEs(landscape, save.plot = 'NLT_logistic-concave-to-convex_fees')


