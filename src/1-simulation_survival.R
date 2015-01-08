if(!exists("prj")){
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project'!")
} else {
  source(file.path(prj, "src", "0-load.R"))
}

#####--------------------------------------------------------------------------
### Simulation 2 - Proportions
### Power
#####------------------------------------
# Settings
# no of simulated datasets
# sample sizes
N <- c(3, 6 ,9)
n_animals <- 10
# proportions in effect groups
pEs <- seq(0.6, 0.95, 0.05)
# both as grid
todo1_p <- expand.grid(N = N, pE = pEs)
  


#####------------------------------------
# simulate data
sims1_p <- NULL
set.seed(1234)
for(i in seq_len(nrow(todo1_p))){
  sims1_p[[i]] <- dosim2(N = todo1_p[i, 'N'], 
                         pC = 0.95, pE = todo1_p[i, 'pE'], 
                         nsim = nsims, n_animals = n_animals)
}

# plot one realisation of simulated data
df <- data.frame(x = sims1_p[[22]]$x, y = sims1_p[[22]]$y[ , 7] / 10)
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
pow_glob_p$variable <- factor(pow_glob_p$variable, unique(pow_glob_p$variable)[c(2,1,3)],
                              c('lm', 'glm', 'np'))

plot_pow_glob_p <- ggplot(pow_glob_p) +
  geom_line(aes(y = value, x = pE, group = variable, linetype = variable)) +
  geom_point(aes(y = value, x = pE, shape = variable), size = 4, color = 'black') +
  facet_grid( ~N, labeller = n_labeller) +  
  # axes
  labs(x = expression('pE'), 
       y = expression(paste('Power (global test , ', alpha, ' = 0.05)'))) +
  # appearance
  mytheme +
  scale_shape_manual('Method', values=c(0,2,4)) +
  scale_linetype_discrete('Method') +
  ylim(c(0,1))
plot_pow_glob_p


### loec
pow_loec_p <- ldply(res1_p, p_loec, type = 'power')
pow_loec_p$pE <- todo1_p$pE
pow_loec_p$N <- todo1_p$N
pow_loec_p  <- melt(pow_loec_p, id.vars = c('pE', 'N'))
pow_loec_p <- pow_loec_p[pow_loec_p$variable %in% c('loeclm', 'loecglm', 'loecpw'), ]
pow_loec_p$variable <- factor(pow_loec_p$variable, unique(pow_loec_p$variable), 
                              labels = c('lm', 'glm', 'np'))

plot_pow_loec_p <- ggplot(pow_loec_p) +
  geom_line(aes(y = value, x = pE, group = variable, linetype = variable)) +
  geom_point(aes(y = value, x = pE, shape = variable), 
            size = 4,  color = 'black') +
  facet_grid( ~N, labeller = n_labeller) + 
  labs(x = 'pE', 
       y = expression(paste('Power (LOEC , ', alpha, ' = 0.05)'))) + 
  mytheme +
  scale_shape_manual('Method', values=c(0,2,4)) +
  scale_linetype_discrete('Method') +
  ylim(c(0,1))
# clean up
# rm(res1_p, df, sims1_p, todo1_p, pEs, n_animals, N)



### ----------------------------------------------------------------------------
# Type1 Error
# Settings
# sample sizes
N <- c(3, 6 ,9)
n_animals <- 10
# proportions 
ps <- seq(0.6, 0.95, 0.05)
# both as grid
todo2_p <- expand.grid(N = N, ps = ps)


#####------------------------------------
# simulate data
sims2_p <- NULL
for(i in seq_len(nrow(todo2_p))){
  sims2_p[[i]] <- dosim2(N = todo2_p[i, 'N'], 
                      pC = todo2_p[i, 'ps'], pE = todo2_p[i, 'ps'],    # pC = CE
                      nsim = nsims, n_animals = n_animals)
}



#####------------------------------------
# analyse data
if(sim2){
  res2_p <- llply(sims2_p, resfoo2, asin = 'ecotox', .progress = 'text')
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
t1_glob_p$variable <- factor(t1_glob_p$variable, unique(t1_glob_p$variable)[c(2,1,3)],
                             c('lm', 'glm', 'np'))


plot_t1_glob_p <- ggplot(t1_glob_p) +
  geom_line(aes(y = value, x = ps, group = variable, linetype = variable)) +
  geom_point(aes(y = value, x = ps, shape = variable), 
             size = 4, color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = 0.05, yend = 0.05), 
               linetype = 'dashed') + 
  facet_grid( ~N, labeller = n_labeller) +  
  # axes
  labs(x = 'p', 
       y = expression(paste('Type 1 error (global test , ', alpha, ' = 0.05)'))) +
  # appearance
  mytheme + 
  # legend title
  scale_shape_manual('Method', values=c(0,2,4)) +
  scale_linetype_discrete('Method') 
plot_t1_glob_p


# LOEC
t1_loec_p <- ldply(res2_p, p_loec, type = 't1')
t1_loec_p$ps <- todo2_p$ps
t1_loec_p$N <- todo2_p$N
t1_loec_p <- melt(t1_loec_p, id.vars = c('ps', 'N'))
t1_loec_p$variable <- factor(t1_loec_p$variable, unique(t1_loec_p$variable), 
                              labels = c('lm', 'glm', 'np'))

plot_t1_loec_p <- ggplot(t1_loec_p) +
  geom_line(aes(y = value, x = ps, group = variable, linetype = variable)) +
  geom_point(aes(y = value, x = ps, shape = variable), 
             size = 4, color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = 0.05, yend = 0.05), 
               linetype = 'dashed') + 
  facet_grid( ~N, labeller = n_labeller) + 
  # axes
  labs(x = 'p', 
       y = expression(paste('Type 1 error (LOEC , ', alpha, ' = 0.05)'))) + 
  # appearance
  mytheme + 
  # legend title
  scale_shape_manual('Method', values=c(0,2,4)) +
  scale_linetype_discrete('Method')
plot_t1_loec_p

# clean up
# rm(res2_p, sims2_p, todo2_p, ps, n_animals, N)
