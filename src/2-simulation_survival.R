#### TODO: check binomial model!
if(!exists('ld')){
  source("/home/edisz/Documents/Uni/Projects/PHD/6USETHEGLM/src/0-load.R")
}
#####--------------------------------------------------------------------------
### Simulation 2 - Proportions
### Power
#####------------------------------------
# Settings
# no of simulated datasets
nsims <- 100
# sample sizes
N <- c(3, 6 ,9 ,12)
n_animals <- 10
# proportions in effect groups
pEs <- seq(0.5, 0.9, 0.05)
# both as grid
todo1_p <- expand.grid(N = N, pE = pEs)
  


#####------------------------------------
# simulate data
sims1_p <- NULL
set.seed(1234)
for(i in seq_len(nrow(todo1_p))){
  sims1_p[[i]] <- dosim2(N = todo1_p[i, 'N'], pE = todo1_p[i, 'pE'], nsim = nsims, n_animals = n_animals)
}

# plot one realisation of simulated data
df <- data.frame(x = sims1_p[[3]]$x, y = sims1_p[[3]]$y[ , 2] / 10)
ggplot(df, aes(x = x, y = y)) +
  geom_boxplot(fill = 'grey80') +
  scale_x_discrete(labels = c('C', 'T1', 'T2', 'T3', 'T4', 'T5')) +
  labs(x = 'Treatment', y = 'Prop. surv.') +
  theme_bw(base_size = 12, 
           base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),
        text = element_text(size=14),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))



#####------------------------------------
# analyse simulations
if(sim2){
  res1_p <- llply(sims1_p, resfoo2, .progress = 'text')
  saveRDS(res1_p, file.path(cachedir, 'res1_p.rds'))
} else {
  res1_p <- readRDS(file.path(cachedir, 'res1_p.rds'))
}


#####------------------------------------
# Results
# global power
pow_glob_p <- ldply(res1_p, p_glob)
pow_glob_p$pE <- todo1_p$pE
pow_glob_p$N <- todo1_p$N
pow_glob_p <- melt(pow_glob_p, id.vars = c('pE', 'N'))

plot_pow_glob_p <- ggplot(pow_glob_p) +
  geom_line(aes(y = value, x = pE, group = variable)) +
  geom_point(aes(y = value, x = pE, fill = variable), size = 4, pch = 21, color = 'black') +
  coord_trans(xtrans = 'log10') +
  # use here
  scale_x_continuous(breaks = pEs) +
  facet_wrap(~N) + 
  # axes
  labs(x = expression('pE'), 
       y = expression(paste('Power (global test , ', alpha, ' = 0.05)'))) +
  # appearance
  theme_bw(base_size = 12, 
           base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        text = element_text(size=14),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  # legend
  scale_fill_grey(name = '', 
                  breaks = c('lm', 'glm', 'k'), 
                  labels = c('LM + arcsin', 'GLM (bin.)', 'Kruskal'),
                  start = 0, end = 1) +
  theme(legend.position="bottom", legend.key = element_blank())
plot_pow_glob_p


### loec
pow_loec_p <- ldply(res1_p, p_loec, type = 'power')
pow_loec_p$pE <- todo1_p$pE
pow_loec_p$N <- todo1_p$N
pow_loec_p <- melt(pow_loec_p, id.vars = c('pE', 'N'))

plot_pow_loec_p <- ggplot(pow_loec_p) +
  geom_line(aes(y = value, x = pE, group = variable)) +
  geom_point(aes(y = value, x = pE, fill = variable), size = 4, pch = 21, color = 'black') +
  facet_grid( ~N) + 
  # axes
  labs(x = 'pE', 
       y = expression(paste('Power (LOEC , ', alpha, ' = 0.05)'))) + 
  # appearance
  theme_bw(base_size = 12, 
           base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        text = element_text(size = 14),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 14,face = "bold")) +
  # legend
  scale_fill_grey(name = '', 
                  breaks = c('lm', 'glm', 'pw'), 
                  labels = c('LM + arcsin', 'GLM (bin.)', 'Wilcox'),
                  start = 0, end = 1) +
  theme(legend.position="bottom", legend.key = element_blank())
plot_pow_loec_p

# clean up
rm(res1_p, df, sims1_p, todo1_p, pEs, n_animals, N, nsims)



### ----------------------------------------------------------------------------
# Type1 Error
# Settings
# no of simulated datasets
nsims <- 100
# sample sizes
N <- c(3, 6 ,9 ,12)
n_animals <- 10
# proportions 
ps <- seq(0.5, 0.9, 0.05)
# both as grid
todo2_p <- expand.grid(N = N, ps = ps)


#####------------------------------------
# simulate data
sims2_p <- NULL
set.seed(1234)
for(i in seq_len(nrow(todo2_p))){
  sims2_p[[i]] <- dosim2(N = todo2_p[i, 'N'], 
                      pC = todo2_p[i, 'ps'], pE = todo2_p[i, 'ps'],    # pC = CE
                      nsim = nsims, n_animals = n_animals)
}


#####------------------------------------
# analyse data
if(sim2){
  res2_p <- llply(sims2_p, resfoo2, .progress = 'text')
  saveRDS(res2_p, file.path(cachedir, 'res2_p.rds'))
} else {
  res2_p <- readRDS(file.path(cachedir, 'res2_p.rds'))
}


#####------------------------------------
# Results
# Global test (how often wrongly assigned an effect)
t1_glob_p <- ldply(res2_p, p_glob)
t1_glob_p$ps <- todo2_p$ps
t1_glob_p$N <- todo2_p$N
t1_glob_p <- melt(t1_glob_p, id.vars = c('ps', 'N'))

plot_t1_glob_p <- ggplot(t1_glob_p) +
  geom_line(aes(y = value, x = ps, group = variable)) +
  geom_point(aes(y = value, x = ps, fill = variable), size = 4, pch = 21, color = 'black') +
  geom_segment(aes(x = .5, xend = 0.99, y = 0.05, yend = 0.05), linetype = 'dashed') + 
  coord_trans(xtrans = 'log2') +
  facet_grid(~N) + 
  # axes
  labs(x = 'pE', 
       y = expression(paste('Type 1 error (global test , ', alpha, ' = 0.05)'))) +
  # appearance
  theme_bw(base_size = 12, 
           base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        text = element_text(size=14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14,face="bold")) +
  # legend
  scale_fill_grey(name = '', 
                  breaks = c('lm', 'glm', 'pk'), 
                  labels = c('LM + arcsin', 'GLM (neg. bin.)', 'Kruskal'),
                  start = 0, end = 1) +
  theme(legend.position = "bottom", legend.key = element_blank())
plot_t1_glob_p


# LOEC
t1_loec_p <- ldply(res2_p, p_loec, type = 't1')
t1_loec_p$ps <- todo2_p$ps
t1_loec_p$N <- todo2_p$N
t1_loec_p <- melt(t1_loec_p, id.vars = c('ps', 'N'))

plot_t1_loec_p <- ggplot(t1_loec_p) +
  geom_line(aes(y = value, x = ps, group = variable)) +
  geom_point(aes(y = value, x = ps, fill = variable), size = 4, pch = 21, color = 'black') +
  geom_segment(aes(x = .5, xend = 0.99, y = 0.05, yend = 0.05), linetype = 'dashed') + 
  coord_trans(xtrans = 'log2') +
  facet_grid(~N) + 
  # axes
  labs(x = 'pE', 
       y = expression(paste('Type 1 error (LOEC , ', alpha, ' = 0.05)'))) + 
  # appearance
  theme_bw(base_size = 12, 
           base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        text = element_text(size = 14),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 14, face = "bold")) +
  # legend
  scale_fill_grey(name = '', 
                  breaks = c('lm', 'glm', 'pw'), 
                  labels = c('LM + arcsin', 'GLM (neg. bin.)', 'Wilcox'),
                  start = 0, end = 1) +
  theme(legend.position="bottom", legend.key = element_blank())
plot_t1_loec_p

# clean up
rm(res2_p, df, sims2_p, todo2_p, ps, n_animals, N, nsims)
