rm(list = ls())

# load auxiliary functions
source('./ecoFunctions.R')
library(tidyr)
library(CEGO)

# tunable parameters
n_mut <- 6
mycolors <- c('#939598', '#be1e2d', '#85c441', '#d68f28', '#415ba9', '#a96cad')
mycolors <- c('#939598', '#d68f28', '#415ba9', '#a96cad')

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
  
  dy <- min(c(max(landscape$f) - landscape$f[1], landscape$f[1] - min(landscape$f)))
  dy <- round(dy/0.1)*0.1
  ybreaks <- seq(landscape$f[1] - 10*dy, landscape$f[1] + 10*dy, by = dy)
  
  myplot <-
    ggplot(df, aes(x = nmut, y = f, group = edge_id, color = color)) +
      geom_abline(slope = 0,
                  intercept = landscape$f[landscape$genot == ''],
                  color = '#d1d3d4') +
      geom_line() +
      scale_x_continuous(name = '# mutations',
                         breaks = 0:n_mut,
                         labels = as.character(0:n_mut)) +
      scale_y_continuous(name = 'Phenotype (a.u.)',
                         breaks = ybreaks,
                         expand = c(0.05, 0.05)) +
      scale_color_manual(values = setNames(mycolors, LETTERS[1:length(mycolors)])) +
      theme_bw() +
      theme(aspect.ratio = 0.6,
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = 'none',
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 16),
            axis.ticks.length = unit(-0.15, 'cm')) +
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)
  
  if (save.plot != F) {
    ggsave(myplot,
           file = paste('../plots/', save.plot, '.pdf', sep = ''),
           dpi = 600,
           width = 100,
           height = 80,
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
  
  dy <- 0.85*min(c(max(gedf$d_f), -min(gedf$d_f)))
  dy <- round(dy/0.1)*0.1
  ybreaks <- seq(-10*dy, 10*dy, by = dy)
  
  dx <- 0.85*(max(gedf$background_f) - min(gedf$background_f))/2
  dx <- round(dx/0.1)*0.1
  xm <- round(mean(gedf$background_f)/0.1)*0.1
  xbreaks <- seq(xm - 10*dx, xm + 10*dx, by = dx)
  
  myplot <-
    ggplot(gedf, aes(x = background_f, y = d_f, group = knock_in, color = color)) +
      geom_abline(slope = 0,
                  intercept = 0,
                  color = '#d1d3d4') +
      geom_point() +
      geom_smooth(method = 'lm',
                  formula = y~x,
                  se = F,
                  fullrange = T,
                  color = 'black') +
      facet_wrap(~knock_in) +
      scale_x_continuous(name = expression(paste(italic(F), ' (background) (a.u.)', sep = '')),
                         breaks = xbreaks,
                         expand = c(0.05, 0.05)) +
      scale_y_continuous(name = expression(paste(Delta, italic(F), ' (a.u.)', sep = '')),
                         breaks = ybreaks,
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
            axis.text = element_text(size = 16),
            axis.ticks.length = unit(-0.15, 'cm'),
            axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)
  
  if (save.plot != F) {
    ggsave(myplot,
           file = paste('../plots/', save.plot, '.pdf', sep = ''),
           dpi = 600,
           width = 140,
           height = 110,
           units = 'mm')
  }
  
  return(myplot)
  
}

# build landscape with only pairwise interactions
makePairwiseLandscape <- function(eps) {
  
  f <- setNames(rep(NA, length(genots)),
                genots)
  f[1] <- 1
  for (g in genots[2:length(genots)]) {
    
    gs <- as.numeric(strsplit(g, split = ',')[[1]])
    f[g] <- 1 + sum(diag(eps)[gs]) # additive contribution
    
    if (length(gs) > 1) { # epistatic contribution
      ps <- t(combn(gs, 2))
      for (i in 1:nrow(ps)) {
        f[g] <- f[g] + eps[ps[i,1], ps[i,2]]
      }
    }
    
  }
  
  return(data.frame(genot = genots,
                    f = as.numeric(f)))

}

### LANDSCAPES

# 0: random landscape
landscape <- data.frame(genot = genots,
                        f = rnorm(length(genots), mean = 1, sd = 0.2))
landscape$f[1] <- 1

plotGraph(landscape, save.plot = 'unf_random_graph')
plotFEEs(landscape, save.plot = 'unf_random_fees')

# 1: additive landscape
set.seed(0)
s_additive <- c(rnorm(3, mean = 0.1, sd = 0.05),
                rnorm(3, mean = -0.1, sd = 0.05))
eps <- matrix(0, nrow = n_mut, ncol = n_mut)
diag(eps) <- s_additive
landscape <- makePairwiseLandscape(eps)

plotGraph(landscape, save.plot = 'unf_additive_graph')
plotFEEs(landscape, save.plot = 'unf_additive_fees')

# 2: positive epistasis between 1-2
eps <- matrix(0, nrow = n_mut, ncol = n_mut)
diag(eps) <- s_additive
eps[1,2] <- 0.1
landscape <- makePairwiseLandscape(eps)
landscape$color <- c('A', 'B')[1 + (grepl('1', landscape$genot) & grepl('2', landscape$genot))]

plotGraph(landscape, save.plot = 'unf_positiveEpist-1-2_graph')
plotFEEs(landscape, save.plot = 'unf_positiveEpist-1-2_fees')

# 3: negative magnitude epistasis between 1-2
eps <- matrix(0, nrow = n_mut, ncol = n_mut)
diag(eps) <- s_additive
eps[1,2] <- -0.25*(s_additive[1] + s_additive[2])
landscape <- makePairwiseLandscape(eps)
landscape$color <- c('A', 'B')[1 + (grepl('1', landscape$genot) & grepl('2', landscape$genot))]

plotGraph(landscape, save.plot = 'unf_negativeMagnEpist-1-2_graph')
plotFEEs(landscape, save.plot = 'unf_negativeMagnEpist-1-2_fees')

# 4: negative sign epistasis between 1-2
eps <- matrix(0, nrow = n_mut, ncol = n_mut)
diag(eps) <- s_additive
eps[1,2] <- -1.5*((s_additive[1] + s_additive[2]) - max(c(s_additive[1], s_additive[2])))
landscape <- makePairwiseLandscape(eps)
landscape$color <- c('A', 'B')[1 + (grepl('1', landscape$genot) & grepl('2', landscape$genot))]

plotGraph(landscape, save.plot = 'unf_negativeSignEpist-1-2_graph')
plotFEEs(landscape, save.plot = 'unf_negativeSignEpist-1-2_fees')

# 5: positive epistasis between 1-2 and between 1-3
eps <- matrix(0, nrow = n_mut, ncol = n_mut)
diag(eps) <- s_additive
eps[1,2] <- 0.1
eps[1,3] <- 0.2
landscape <- makePairwiseLandscape(eps)
landscape$color <- 'A'
landscape$color[grepl('1', landscape$genot) & grepl('2', landscape$genot)] <- 'B'
landscape$color[grepl('1', landscape$genot) & grepl('3', landscape$genot)] <- 'C'
landscape$color[grepl('1', landscape$genot) & grepl('2', landscape$genot) & grepl('3', landscape$genot)] <- 'D'

plotGraph(landscape, save.plot = 'unf_positiveEpist-1-2-and-1-3_graph')
plotFEEs(landscape, save.plot = 'unf_positiveEpist-1-2-and-1-3_fees')

# 6: positive epistasis between 1-2, negative magnitude epistasis between 1-3
eps <- matrix(0, nrow = n_mut, ncol = n_mut)
diag(eps) <- s_additive
eps[1,2] <- -0.05
eps[1,3] <- 0.08
landscape <- makePairwiseLandscape(eps)
landscape$color <- 'A'
landscape$color[grepl('1', landscape$genot) & grepl('2', landscape$genot)] <- 'B'
landscape$color[grepl('1', landscape$genot) & grepl('3', landscape$genot)] <- 'C'
landscape$color[grepl('1', landscape$genot) & grepl('2', landscape$genot) & grepl('3', landscape$genot)] <- 'D'

plotGraph(landscape, save.plot = 'unf_positiveEpist-1-2-negativeMagnEpist-1-3_graph')
plotFEEs(landscape, save.plot = 'unf_positiveEpist-1-2-negativeMagnEpist-1-3_fees')

# 7: positive epistasis between 1-6
eps <- matrix(0, nrow = n_mut, ncol = n_mut)
diag(eps) <- s_additive
eps[1,6] <- 0.1
landscape <- makePairwiseLandscape(eps)
landscape$color <- 'A'
landscape$color[grepl('1', landscape$genot) & grepl('6', landscape$genot)] <- 'B'

plotGraph(landscape, save.plot = 'unf_positiveEpist-1-6_graph')
plotFEEs(landscape, save.plot = 'unf_positiveEpist-1-6_fees')

# 8: megative epistasis between 1-6
eps <- matrix(0, nrow = n_mut, ncol = n_mut)
diag(eps) <- s_additive
eps[1,6] <- -0.1
landscape <- makePairwiseLandscape(eps)
landscape$color <- 'A'
landscape$color[grepl('1', landscape$genot) & grepl('6', landscape$genot)] <- 'B'

plotGraph(landscape, save.plot = 'unf_negativeEpist-1-6_graph')
plotFEEs(landscape, save.plot = 'unf_negativeEpist-1-6_fees')

