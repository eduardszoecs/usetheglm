if(!exists("prj")){
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project'!")
} else {
  source(file.path(prj, "src", "0-load.R"))
}

#####--------------------------------------------------------------------------
##### Simulation 1 -  Count data
### Power
#####------------------------------------
# Settings
# no of simulated datasets
nsims <- 100
# sample sizes
N <- c(3, 6 ,9 ,12)
ctrl <- 2^(c(1:5, 7, 9))
# both as grid
todo1_c <- expand.grid(N = N, ctrl = ctrl)
# fixed theta
theta  <- rep(3.91, 6)  


#####------------------------------------
# simulate data
sims1_c <- NULL
set.seed(seed)
for(i in seq_len(nrow(todo1_c))){
  N <- todo1_c[i, 'N']
  takectrl <- todo1_c[i, 'ctrl']
  # reduce t2-t5 to 50%
  taketrt <- takectrl * 0.5
  mu <- c(rep(takectrl, each = 2), rep(taketrt, each = 4))
  sims1_c[[i]] <- dosim1(N = N, mu = mu, nsims = nsims, theta = theta)
}

# plot one realisation of simulated data
todo1_c[15, ]
df <- data.frame(x = sims1_c[[15]]$x, y = sims1_c[[15]]$y[ , 2])
df$yt <- log(1 / min(df$y[df$y != 0]) * df$y + 1)
dfm <- melt(df)
levels(dfm$variable) <- c('y', 'ln(Ay + 1)')
ggplot(dfm, aes(x = x, y = value)) +
  geom_boxplot(fill = 'grey80') +
  facet_wrap( ~variable, scales = 'free_y') +
  scale_x_discrete(labels = c('C', 'T1', 'T2', 'T3', 'T4', 'T5')) +
  labs(x = 'Treatment', y = 'Abundance') +
  theme_bw(base_size = 12, 
           base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),
        text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14,face = "bold"))



#####------------------------------------
# analyse simulations
if(sim1){
  res1_c <- llply(sims1_c, resfoo1, .progress = 'text')
  saveRDS(res1_c, file.path(cachedir, 'res1_c.rds'))
} else {
  res1_c <- readRDS(file.path(cachedir, 'res1_c.rds'))
}


#####------------------------------------
# Results
# global power
pow_glob_c <- ldply(res1_c, p_glob1)
pow_glob_c$muc <- todo1_c$ctrl
pow_glob_c$N <- todo1_c$N
# restructure data
pow_glob_c <- reshape(pow_glob_c, direction = 'long', idvar = 7:8,
        varying = list(1:3, 4:6),
        v.names = c('power', 'conv'),
        timevar = 'variable',
        times = c('lm', 'glm', 'pk'))

plot_pow_glob_c <- ggplot(pow_glob_c) +
  geom_line(aes(y = power, x = muc, group = variable)) +
  geom_point(aes(y = power, x = muc, fill = variable, size = conv), 
             pch = 21, color = 'black') +
  coord_trans(xtrans = 'log2') +
  scale_x_continuous(breaks = round(unique(todo1_c$ctrl), 0)) +
  facet_grid( ~ N) + 
  # axes
  labs(x = expression(mu[C]), 
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
  scale_fill_grey(name = 'Method', 
                  breaks = c('lm', 'glm', 'pk'), 
                  labels = c('LM + log(Ay+1)', 'GLM (neg. bin.)', 'Kruskal'),
                  start = 0, end = 1) +
  scale_size_continuous( 'Convergence',
    range = c(1,4),
    breaks = seq(70, 100, 10),
    labels = paste0(seq(70, 100, 10), ' %')
    ) +
  theme(legend.position="bottom")
plot_pow_glob_c


### loec
pow_loec_c <- ldply(res1_c, p_loec1, type = 'power')
pow_loec_c$muc <- todo1_c$ctrl
pow_loec_c$N <- todo1_c$N
pow_loec_c  <- reshape(pow_loec_c, direction = 'long', idvar = 7:8,
                       varying = list(1:3, 4:6),
                       v.names = c('power', 'conv'),
                       timevar = 'variable',
                       times = c('lm', 'glm', 'pk'))

plot_pow_loec_c <- ggplot(pow_loec_c) +
  geom_line(aes(y = power, x = muc, group = variable)) +
  geom_point(aes(y = power, x = muc, fill = variable, size = conv), 
             pch = 21, color = 'black') +
  coord_trans(xtrans = 'log2') +
  scale_x_continuous(breaks = round(unique(todo1_c$ctrl), 0)) +
  facet_grid( ~N) + 
  # axes
  labs(x = expression(mu[C]), 
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
  scale_fill_grey(name = 'Method', 
                  breaks = c('lm', 'glm', 'pk'), 
                  labels = c('LM + log(Ay+1)', 'GLM (neg. bin.)', 
                             'Wilcox'),
                  start = 0, end = 1) +
  scale_size_continuous( 'Convergence',
                         range = c(1,4),
                         breaks = seq(70, 100, 10),
                         labels = paste0(seq(70, 100, 10), ' %')
  ) +
  theme(legend.position="bottom")
plot_pow_loec_c

# clean up
rm(res1_c, todo1_c, df, dfm, sims1_c, takectrl, taketrt, theta, 
   mu, ctrl, nsims, N)



### ----------------------------------------------------------------------------
# Type1 Error
#####------------------------------------
# settings
nsims <- 100
# sample sizes
N <- c(3, 6, 9, 12)
ctrl <- 2^(c(1:5, 7, 9))
# both as grid
todo2_c <- expand.grid(N = N, ctrl = ctrl)
theta  <- rep(3.91, 6)  


#####------------------------------------
# simulate data
sims2_c <- NULL
set.seed(seed)
for(i in seq_len(nrow(todo2_c))){
  N <- todo2_c[i, 'N']
  takectrl <- todo2_c[i, 'ctrl']
  # all treatments with same mean
  taketrt <- takectrl * 1
  mu <- c(rep(takectrl, each = 2), rep(taketrt, each = 4))
  sims2_c[[i]] <- dosim1(N = N, mu = mu, nsims = nsims, theta = theta)
}


#####------------------------------------
# analyse data
if(sim1){
  res2_c <- llply(sims2_c, resfoo1, .progress = 'text')
  saveRDS(res2_c, file.path(cachedir, 'res2_c.rds'))
} else {
  res2_c <- readRDS(file.path(cachedir, 'res2_c.rds'))
}


#####------------------------------------
# Results
# Global test (how often wrongly assigned an effect)
t1_glob_c <- ldply(res2_c, p_glob1)
t1_glob_c$muc <- todo2_c$ctrl
t1_glob_c$N <- todo2_c$N
# restructure data
t1_glob_c <- reshape(t1_glob_c, direction = 'long', idvar = 7:8,
                      varying = list(1:3, 4:6),
                      v.names = c('t1', 'conv'),
                      timevar = 'variable',
                      times = c('lm', 'glm', 'pk'))

plot_t1_glob_c <- ggplot(t1_glob_c) +
  geom_line(aes(y = t1, x = muc, group = variable)) +
  geom_point(aes(y = t1, x = muc, fill = variable, size = conv), 
             pch = 21, color = 'black') +
  geom_segment(aes(x = 2, xend = 1024, y = 0.05, yend = 0.05), 
               linetype = 'dashed') + 
  # maby try log2 transformation?
  coord_trans(xtrans = 'log2') +
  # use here
  scale_x_continuous(breaks = round(ctrl, 0)) +
  facet_grid(~N) + 
  # axes
  labs(x = expression(mu[C]), 
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
                  labels = c('LM + log(Ay+1)', 'GLM (neg. bin.)', 'Kruskal'),
                  start = 0, end = 1) +
  theme(legend.position = "bottom", legend.key = element_blank())  +
  scale_size_continuous( 'Convergence',
                         range = c(1,4),
                         breaks = seq(70, 100, 10),
                         labels = paste0(seq(70, 100, 10), ' %')
  )
plot_t1_glob_c


# LOEC
t1_loec_c <- ldply(res2_c, p_loec1, type = 't1')
t1_loec_c$muc <- todo2_c$ctrl
t1_loec_c$N <- todo2_c$N
t1_loec_c <- reshape(t1_loec_c, direction = 'long', idvar = 7:8,
                     varying = list(1:3, 4:6),
                     v.names = c('t1', 'conv'),
                     timevar = 'variable',
                     times = c('lm', 'glm', 'pk'))

plot_t1_loec_c <- ggplot(t1_loec_c) +
  geom_line(aes(y = t1, x = muc, group = variable)) +
  geom_point(aes(y = t1, x = muc, fill = variable, size = conv), 
             pch = 21, color = 'black') +
  geom_segment(aes(x = 2, xend = 1024, y = 0.05, yend = 0.05), 
               linetype = 'dashed') + 
  coord_trans(xtrans = 'log2') +
  scale_x_continuous(breaks = round(unique(todo2_c$ctrl), 0)) +
  facet_grid(~N) + 
  # axes
  labs(x = expression(mu[C]), 
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
                  labels = c('LM + log(Ay+1)', 'GLM (neg. bin.)', 'Wilcox'),
                  start = 0, end = 1) +
  theme(legend.position="bottom", legend.key = element_blank()) +
  scale_size_continuous( 'Convergence',
                         range = c(1,4),
                         breaks = seq(70, 100, 10),
                         labels = paste0(seq(70, 100, 10), ' %')
  )
plot_t1_loec_c


# clean up
rm(res2_c, todo2_c, sims2_c, takectrl, taketrt, theta, 
   mu, ctrl, nsims, N)
