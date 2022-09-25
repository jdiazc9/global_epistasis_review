library(readr)
library(gtools)
library(tidyverse)

#source("~/global_epistasis_clean/Code/general_functions.R")
#not sourcing, instead just pulling two functions from this

###################
#### functions ####
###################
rowcol<-function(ix,n) { #given index, return row and column
  nr=ceiling(n-(1+sqrt(1+4*(n^2-n-2*ix)))/2)
  nc=n-(2*n-nr+1)*nr/2+ix+nr
  cbind(nr,nc)
}

get_effect_df_dist <- function(df){
  
  communities <- apply(df %>% dplyr::select(-fitness), 1, function(row) paste0(row, collapse = '-')) 
  df <- df %>% mutate(community = communities)
  
  #need to average out duplicates before creating effect_df 
  df <- df %>% group_by(community) %>% mutate(fitness = mean(fitness)) %>% ungroup() %>% distinct(community, fitness, .keep_all = TRUE)
  
  ###
  E <- df %>% dplyr::select(-c(fitness,community)) 
  all_species <- colnames(E) 
  all_lens <- rowSums(E)
  
  ###
  all_dist <- dist(E)
  n_exp <- nrow(E)
  loo <- which(all_dist == 1)
  
  effect_df <- tibble()
  for (i in 1:length(loo)){
    rc <-rowcol(loo[i], n_exp)
    row <- rc[1]
    col <- rc[2]
    
    species <- which(E[row,] - E[col,] != 0)
    if ((E[row,] - E[col,])[species] == 1){
      new_com <- row
      old_com <- col
    } else {
      new_com <- col
      old_com <- row
    }
    tmp <- list(species = all_species[species], new_fitness = df[new_com,]$fitness, 
                baseline_fitness = df[old_com,]$fitness, new_community = paste0(E[new_com,], collapse = '-'), 
                baseline_community = paste0(E[old_com,], collapse = '-'))
    effect_df <- rbind(effect_df, tmp)
  }
  
  if (!any(all_lens == 0)){
    #need to add monocultures, assume f(0) = 0 here
    #if wildtype/ 0 community in df, this will already be taken care of 
    monocultures <- which(all_lens == 1)
    for (i in monocultures){
      species <- which(E[i,] == 1)
      tmp <- list(species = all_species[species], new_fitness = df[i,]$fitness, 
                  baseline_fitness = 0, new_community = paste0(E[i,], collapse = '-'), 
                  baseline_community = paste(rep('0', length(all_species)), collapse = '-'))
      effect_df <- rbind(effect_df, tmp)
    }
  }
  
  effect_df <- effect_df %>% group_by(species) %>% mutate(delta_fitness = new_fitness - baseline_fitness)
  
  return(effect_df)
}

get_linear_fits <- function(effect_df){
  effect_df <- effect_df %>% group_by(species) %>% 
    mutate(intercept = coef(lm(delta_fitness ~ baseline_fitness))[1], slope = coef(lm(delta_fitness ~baseline_fitness))[2])
  
  return(effect_df)
}

########################
###### Read data #######
########################

folder <- c("") #Folder where input file is 
data.file <- c("genotype_fitness.txt") #input text file - no spaces or missing values - tab delimited - 1 header row
open.file <- paste(folder, data.file, sep="")
df <- read.table(open.file, header=T)

##extract genetic data
df$genotype=as.character(df$genotype)
all_genes= sort(c("r", "t", "s", "g", "p"))
N <- length(all_genes)

#rename fitness column 
cols <- colnames(df)
cols[4] <- 'fitness'
colnames(df) <- cols

#pull out community composition and diversity variables
#rename ancestral strain for convenience
df$genotype[which(df$genotype == 'anc')] <- ''
communities <- strsplit(df$genotype, split = '')

#Collect data into nicer form and plot FEEs
#make presence absence
presence <- matrix(0, nrow = nrow(df), ncol = N)
for (i in 1:nrow(df)){
  presence[i, which(all_genes %in% communities[[i]])] <- 1
}

df <- data.frame(presence, fitness = df$fitness)
colnames(df) <- c(all_genes, 'fitness')

#avg out replicates
df <- df %>% group_by(across(c(-fitness))) %>% mutate(fitness = mean(fitness)) %>% distinct(fitness, .keep_all = TRUE) %>% ungroup()

effect_df <- get_effect_df_dist(df)

#plot FEEs
p <- effect_df %>% ggplot(aes(x = baseline_fitness, y = delta_fitness)) + geom_point() + 
  facet_wrap(~species) + geom_hline(yintercept = 0) + theme_bw() + xlab("Baseline Fitness") + ylab("Delta Fitness")
show(p)

#make presence absence (convenient later)
presence <- matrix(0, nrow = nrow(effect_df), ncol = N)
coms <- strsplit(effect_df$baseline_community, '-')
for (i in 1:nrow(effect_df)){
  presence[i, ] <- as.numeric(coms[[i]])
}

colnames(presence) <- colnames(df)[1:N]

effect_df <- data.frame(effect_df, presence)

# look at fee_fits
fee_fits <- get_linear_fits(effect_df) %>% distinct(species, intercept, slope)
fee_fits

#################################
## get epistatic coefficients ###
#################################

#combinatorially complete so we can just explicitly solve for these 
y <- df$fitness
f <- as.formula(y ~ .^5)
x <- model.matrix(f, df %>% select(-fitness))

epi_coefs <- MASS::ginv(x) %*% df$fitness

#arrange
epi_coefs <- data.frame(effect = colnames(x), value = epi_coefs)

epi_presence <- matrix(0, nrow = nrow(epi_coefs), ncol = N)
genotypes <- strsplit(epi_coefs$effect, ":")

for (i in 1:length(genotypes)){
  epi_presence[i, which(all_genes %in% genotypes[[i]])] <- 1
}

order <- rowSums(epi_presence) #order of interaction
epi_coefs <- data.frame(epi_presence,epi_coefs, order)
colnames(epi_coefs) <- c(all_genes, 'effect', 'value', 'order')

#plot distribution of coefs by order
#zeroth order is intercept ,first order = additive effect, second order = /epsilon_ij, etc
epi_means <- epi_coefs %>% group_by(order) %>% summarize(mean_value = mean(value))
epi_coefs %>% ggplot(aes(x = value)) + geom_histogram(binwidth = .05) + geom_vline(data = epi_means, aes(xintercept = mean_value)) + facet_wrap(~order) + theme_bw() + ggtitle("Distribution of epistatic coefficients by order")

#The vertical line is the mean of these coefficients, and for order > 1 (ie $\epsilon_{ij}, \epsilon_{ijk}$, etc) 
#the mean seems to be around 0. 

########################################################
#### Doing the same for the fourier coefficients... ####
########################################################

#create new basis 
fourier_basis <- df %>% select(-fitness)
fourier_basis[fourier_basis == 0] <- -1 

fourier_df <- data.frame(fourier_basis, fitness = df$fitness) 

#solve for fourier coefficients
y <- fourier_df$fitness
f <- as.formula(y ~ .^5)
x <- model.matrix(f, fourier_df %>% select(-fitness))

fourier_coefs <- MASS::ginv(x) %*% df$fitness

#collect into convenient form 
fourier_coefs <- data.frame(effect = colnames(x), value = fourier_coefs)

fourier_presence <- matrix(0, nrow = nrow(fourier_coefs), ncol = N)
genotypes <- strsplit(fourier_coefs$effect, ":")

for (i in 1:length(genotypes)){
  fourier_presence[i, which(all_genes %in% genotypes[[i]])] <- 1
}

#process for plotting
order <- rowSums(fourier_presence)
fourier_coefs <- data.frame(fourier_presence,fourier_coefs, order )
colnames(fourier_coefs) <- c( all_genes,'effect', 'value', 'order')

#plot distribution of coefficients
fourier_means <- fourier_coefs %>% filter(order != 0) %>% group_by(order) %>% summarize(mean_value = mean(value))
fourier_coefs  %>% filter(order != 0) %>% ggplot(aes(x = value)) + geom_histogram() + geom_vline(data = fourier_means, aes(xintercept = mean_value)) + facet_wrap(~order) + theme_bw() + ggtitle("Distribution of fourier coefficients by order")


#We again see the distribution centered around zero for orders > 1

##########################
### get v_i ##############
##########################
#We can get $\tilde{v}_i$ and thus the predicted slope/intercept of the FEEs from the fourier coefficients

#ugly loop but not sure a nested list is really better... 
vi <- c()
for (i in 1:N){
  focal_mut <- colnames(df)[i]
  focal_fi <- fourier_coefs$value[which(fourier_coefs$effect == focal_mut)] #get individual effect
  focal_fijs <- fourier_coefs$value[which(fourier_coefs[,i] == 1 & fourier_coefs$order == 2)] 
  
  denom <- 0
  for (j in (1:N)[-i]){
    mut_j <- colnames(df)[j]
    fj <- fourier_coefs$value[which(fourier_coefs$effect == mut_j)]
    mut_j_fij <- fourier_coefs$value[which(fourier_coefs[,i] == 1 & fourier_coefs[,j] == 1 & fourier_coefs$order == 2)]
    denom <- denom + (fj - mut_j_fij)^2
  }
  vi[i] <- (sum(focal_fijs^2) - sum(focal_fi * focal_fijs))/(denom)
}

#collect
vis <- data.frame(species = colnames(df)[1:N], vi = vi)

##############################################
##### plot predicted slope vs actual ########
#############################################
to_plot <- inner_join(vis, fee_fits, by = 'species')
p <- ggplot(to_plot, aes(x = slope, y = -2*vi)) + geom_point() + geom_abline() + theme_bw() + ggtitle("Observed FEE slopes vs Predicted slopes") + xlab('Fitted slope')
show(p)


#############################################
### plot avg epsilon ijs and deltaj [-i] ###
#############################################

#the list of epsilon i,j averaged over every background where we can place the i,j pair -- 1/4 fij
#the list of delta yj {-i} the average fitness of mutation j across all backgrounds where i is not present

#ugly loop!
res <- tibble()
for (i in 1:N){
  mut_i <- colnames(df)[i] #name 
  for (j in (1:N)[-i]){
    mut_j <- colnames(df)[j] #name 
    
    #get values 
    focal_Eij <- (4)*fourier_coefs$value[which(fourier_coefs$order == 2 & fourier_coefs[,i] == 1 & fourier_coefs[,j] == 1)]
    focal_delta_j <- effect_df %>% filter(species == colnames(df)[j], get(mut_i) == 0) %>% summarize(mean(delta_fitness))
    focal_fj <- fourier_coefs$value[which(fourier_coefs$order == 1 & fourier_coefs[,j] == 1)]
    
    #append
    tmp <- list(mut_i = mut_i, mut_j = mut_j, Eij = focal_Eij, fj = focal_fj, delta_j = focal_delta_j[[1]], 
                slope = fee_fits$slope[which(fee_fits$species == mut_i)])
    res <- rbind(res, tmp)
  }
}

#add dummy variable so this shows up in plot label
res <- res %>% mutate(mut_i_label = paste0(mut_i,', slope = ', round(slope, 3)), sep = '')
#plot
ggplot(res, aes(x = Eij, y = delta_j, color = mut_j)) + geom_point() + theme_bw() + 
  ggtitle("") + facet_wrap(~mut_i_label) + labs(color = 'Mutation j') + geom_vline(xintercept = 0, alpha = .25) + 
  geom_hline(yintercept = 0, alpha = .25) + xlab('E_ijs averaged over all backgrounds') + ylab('Delta Y_{j}[-i]')

###############################
### moving on to intercepts ###
################################

y_bar <- fourier_coefs$value[1] #fourier intercept 
alpha <- c()
for (i in 1:N){
  f_i <- fourier_coefs$value[which(fourier_coefs[,i] == 1 & fourier_coefs$order == 1)]
  alpha[i] <- 2* f_i* (1 - vi[i]) + 2 * vi[i] * y_bar #def 
}

#collect
alpha_is <- data.frame(species = colnames(df)[1:N], alpha_i = alpha)

alphas_to_plot <- inner_join(alpha_is, fee_fits, by = 'species')

#plot actual vs predicted intercepts
p <- ggplot(alphas_to_plot, aes(x = intercept, y = alpha_i)) + geom_point() + geom_abline() + theme_bw() + ggtitle("Observed FEE intercepts vs predicted intercepts") + xlab('Fitted intercept from FEEs')
show(p)

###################################
# plot predicted linear patterns #
####################################
predicted_val_df <- cbind(vis, alpha)

effect_df <- left_join(effect_df, predicted_val_df, by = 'species') 

#plot just the linear patterns we're predicting 
effect_df %>% ggplot(aes(x = baseline_fitness, y = delta_fitness)) + geom_point() + 
  geom_hline(yintercept = 0) + geom_abline(aes(intercept = alpha, slope = -2*vi), color ='purple') + 
  facet_wrap(~species) + theme_bw() + xlab("Baseline Fitness") + ylab("Delta Fitness") + 
  ggtitle('Effects of mutations with predicted linear pattern overlaid')


#plot the linear patterns with the line of best fit as well  
effect_df %>% ggplot(aes(x = baseline_fitness, y = delta_fitness)) + geom_point() + 
  geom_hline(yintercept = 0) + geom_abline(aes(intercept = alpha, slope = -2*vi), color ='purple') + 
  facet_wrap(~species) + theme_bw() + xlab("Baseline Fitness") + ylab("Delta Fitness") + 
  geom_smooth(method = 'lm') + ggtitle('Effects of mutations with predicted linear pattern overlaid + line of best fit')


########################
# juan's plots
########################

# eps_ij vs deltaF_j

mut_names <- setNames(c('rbs', 'topA', 'spoT', 'glmUS', 'pykF'),
                      c('r', 't', 's', 'g', 'p'))
res$mut_i <- mut_names[res$mut_i]
res$mut_j <- mut_names[res$mut_j]

res$mut_i <- factor(res$mut_i, levels = as.character(mut_names))
res$mut_j <- factor(res$mut_j, levels = as.character(mut_names))

res$prod <- abs(res$Eij*res$delta_j)
res$color <- c('A', 'B')[1 + (res$Eij>0)]

myplot <- 
  ggplot(res, aes(x = Eij, y = delta_j, size = prod, color = color)) +
    geom_vline(xintercept = 0, color = '#d1d3d4') + 
    geom_hline(yintercept = 0, color = '#d1d3d4') +
    geom_point() +
    geom_text(aes(label = mut_j),
              size = 4,
              hjust = -0.5,
              fontface = 'italic') +
    facet_wrap(~mut_i, nrow = 1) +
    scale_color_manual(values = c('#be1e2d', '#85c441')) +
    scale_x_continuous(name = expression(paste('<', italic(epsilon[ij]), '>', sep = '')),
                       expand = c(0.01, 0.01),
                       breaks = c(-0.025, 0, 0.025),
                       labels = c('-0.025', '0', '0.025')) +
    scale_y_continuous(name = expression(paste('<', italic(delta), 'F', '>', sep = '')),
                       expand = c(0.015, 0.015)) +
    theme_bw() +
    theme(aspect.ratio = 0.6,
          panel.grid = element_blank(),
          legend.position = 'none',
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          axis.ticks.length = unit(-0.15, 'cm'),
          strip.background = element_blank(),
          strip.text = element_text(size = 16, face = 'italic'))

print(myplot)

ggsave(myplot,
       file = paste('../plots/', 'fig3_a', '.pdf', sep = ''),
       dpi = 600,
       width = 300,
       height = 80,
       units = 'mm',
       limitsize = F)

# FEEs

effect_df$species <- mut_names[effect_df$species]
effect_df$species <- factor(effect_df$species, levels = as.character(mut_names))

myplot <-
ggplot(effect_df, aes(x = baseline_fitness, y = delta_fitness)) +
  geom_hline(yintercept = 0,
             color = '#d1d3d4') +
  geom_smooth(method = 'lm',
              formula = y~x,
              fullrange = T,
              size = 1,
              color = 'gray') +
  geom_point() +
  geom_abline(aes(intercept = alpha, slope = -2*vi),
              color = 'black',
              size = 1) + 
  facet_wrap(~species, nrow = 1) +
  scale_x_continuous(name = expression(paste(italic(F), ' (background)', sep = '')),
                     breaks = c(1, 1.2)) +
  scale_y_continuous(name = expression(paste(Delta, italic(F), sep = '')),
                     breaks = c(0, 0.1, 0.2)) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        legend.position = 'none',
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.ticks.length = unit(-0.15, 'cm'),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, face = 'italic'))

print(myplot)

ggsave(myplot,
       file = paste('../plots/', 'fig3_b', '.pdf', sep = ''),
       dpi = 600,
       width = 295,
       height = 80,
       units = 'mm',
       limitsize = F)


